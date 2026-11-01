#!/usr/bin/env python3
"""登山道トレース → KML・経路図・標高図

OSM（api.openstreetmap.org）と国土地理院（cyberjapandata.gsi.go.jp）の登山道から
開始点→終了点の最短ルートを引き、KML・経路図PNG・標高図PNG・レポートJSONを出力する。

使い方:
  python3 trailkml.py --mountain 鳥海山 --start 鉾立 --goal 山頂 --out chokai
  python3 trailkml.py --mountain 鳥海山 --start "39°07'18.5\"N 139°59'29.7\"E" --goal 山頂 --out chokai
  地点は「緯度,経度」「度分秒」「地点名（地理院の地名検索→OSMの地点名の順に探す）」「山頂」のいずれか。
  山名は地理院の地名検索で位置を決める（同名の山は標高が最も高いものを使う）。
"""
import argparse, heapq, io, json, math, os, re, shlex, statistics, subprocess, sys, time
VERSION = "2.3"
from concurrent.futures import ThreadPoolExecutor
import requests
from PIL import Image, ImageDraw, ImageFont

UA = {"User-Agent": "mountain-trail-kml/1.0 (personal hiking route tool)"}
GSI = "https://cyberjapandata.gsi.go.jp/xyz"
FONT_DEFAULT = "/usr/share/fonts/opentype/noto/NotoSansCJK-Medium.ttc"
FONT_CANDIDATES = [
    FONT_DEFAULT,
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/noto-cjk/NotoSansCJK-Medium.ttc",
    "/usr/share/fonts/google-noto-cjk/NotoSansCJK-Medium.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W5.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W6.ttc",
    "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
    "/Library/Fonts/NotoSansCJKjp-Medium.otf",
    "/System/Library/Fonts/Hiragino Sans GB.ttc",
]

def find_font():
    """日本語フォントを探す。既定の Linux パス → 代表的なパス → fc-match の順。見つからなければ None"""
    for p in FONT_CANDIDATES:
        if os.path.exists(p): return p
    for pat in ("Noto Sans CJK JP:style=Medium", "Hiragino Sans", ":lang=ja"):
        try:
            r = subprocess.run(["fc-match", "-f", "%{file}", pat], capture_output=True, text=True, timeout=10)
            f = r.stdout.strip()
            if f and os.path.exists(f) and re.search(r"(?i)cjk|hiragino|ヒラギノ|gothic|ipa|takao|jp", f): return f
        except Exception:
            pass
    return None

FONT = find_font()

def font(size):
    if FONT: return ImageFont.truetype(FONT, size)
    try: return ImageFont.load_default(size=size)
    except TypeError: return ImageFont.load_default()
S = requests.Session(); S.headers.update(UA)

# ---------------------------------------------------------------- 基本
def hav(a, b):
    r = math.pi / 180; dl = (b[1]-a[1])*r; dn = (b[0]-a[0])*r
    h = math.sin(dl/2)**2 + math.cos(a[1]*r)*math.cos(b[1]*r)*math.sin(dn/2)**2
    return 12742000*math.asin(math.sqrt(h))

def plen(p): return sum(hav(p[i-1], p[i]) for i in range(1, len(p)))

def cum(p):
    c = [0.0]
    for i in range(1, len(p)): c.append(c[-1] + hav(p[i-1], p[i]))
    return c

def seg_dist(p, a, b):
    lat0 = math.radians(p[1]); k = 111320
    ax, ay = (a[0]-p[0])*k*math.cos(lat0), (a[1]-p[1])*k
    bx, by = (b[0]-p[0])*k*math.cos(lat0), (b[1]-p[1])*k
    dx, dy = bx-ax, by-ay; L = dx*dx + dy*dy
    t = 0 if L == 0 else max(0, min(1, -(ax*dx + ay*dy)/L))
    return math.hypot(ax+t*dx, ay+t*dy)

def dist_to_line(p, line): return min(seg_dist(p, line[i-1], line[i]) for i in range(1, len(line)))

def tile_xy(lon, lat, z):
    n = 2**z; r = math.radians(lat)
    return (lon+180)/360*n, (1-math.log(math.tan(r)+1/math.cos(r))/math.pi)/2*n

def tile_ll(x, y, z):
    n = 2**z
    return x/n*360-180, math.degrees(math.atan(math.sinh(math.pi*(1-2*y/n))))

def parse_point(s):
    """'39.12,139.99' / '39°07'18.5"N 139°59'29.7"E' → (lon, lat)。座標でなければ None"""
    s = s.strip()
    m = re.findall(r"(\d+(?:\.\d+)?)\s*[°度]\s*(\d+(?:\.\d+)?)\s*['′分]\s*(\d+(?:\.\d+)?)?\s*(?:\"|″|”|秒)?\s*([NSEW北南東西])?", s)
    if len(m) == 2:
        vals = []
        for d, mi, se, h in m:
            v = float(d) + float(mi)/60 + float(se or 0)/3600
            if h in ("S", "W", "南", "西"): v = -v
            vals.append((v, h))
        (a, ha), (b, hb) = vals
        if ha in ("E", "W", "東", "西"): a, b = b, a
        return (b, a)
    m = re.match(r"^\s*(-?\d+\.\d+)\s*[, ]\s*(-?\d+\.\d+)\s*$", s)
    if m:
        a, b = float(m.group(1)), float(m.group(2))
        return (b, a) if abs(a) < 90 and abs(b) > 90 or a < b else (a, b)  # 日本では 緯度<経度
    return None

def dijkstra(adj, s, g):
    dist = {s: 0.0}; prev = {}; h = [(0.0, s)]
    while h:
        d, u = heapq.heappop(h)
        if d > dist.get(u, 1e18): continue
        if u == g: break
        for v, w, wid in adj.get(u, ()):
            nd = d + w
            if nd < dist.get(v, 1e18): dist[v] = nd; prev[v] = (u, wid); heapq.heappush(h, (nd, v))
    if g != s and g not in prev: return None
    path, ways, c = [g], [], g
    while c != s: u, wid = prev[c]; ways.append(wid); c = u; path.append(c)
    return path[::-1], ways[::-1]

def components(adj):
    comp = {}
    for s0 in adj:
        if s0 in comp: continue
        comp[s0] = s0; st = [s0]
        while st:
            u = st.pop()
            for v, _, _ in adj.get(u, ()):
                if v not in comp: comp[v] = s0; st.append(v)
    return comp

def snap_pair(adj, pos, A, B):
    """A・B を、互いにつながっている道の上の最も近い点に吸着させる（孤立した駐車場の通路などを避ける）"""
    comp = components(adj)
    def nearest(p, cid=None):
        best = min((i for i in adj if cid is None or comp[i] == cid), key=lambda i: hav(p, pos(i)))
        return best, hav(p, pos(best))
    opts = []
    for first in ("A", "B"):
        if first == "A": a, da = nearest(A); b, db = nearest(B, comp[a])
        else: b, db = nearest(B); a, da = nearest(A, comp[b])
        opts.append((da + db, a, da, b, db))
    _, a, da, b, db = min(opts)
    return a, da, b, db

# ---------------------------------------------------------------- 地名検索（国土地理院）
def gsi_search(q):
    """地理院の地名検索。[(title, lon, lat)]"""
    try:
        r = S.get("https://msearch.gsi.go.jp/address-search/AddressSearch", params={"q": q}, timeout=30)
        return [(f["properties"].get("title", ""), *f["geometry"]["coordinates"]) for f in r.json()]
    except Exception:
        return []

def locate_mountain(name):
    """山名 → (lon, lat, 候補一覧)。同名の山が複数あるときは、標高が最も高いものを選ぶ"""
    allh = gsi_search(name)
    hits = [h for h in allh if h[0] == name] or [h for h in allh if name in h[0]]
    if not hits: return None
    cands = []
    for t, lo, la in hits[:10]:
        z = elev(lo, la)
        cands.append({"title": t, "lon": lo, "lat": la, "ele": z or 0})
    cands.sort(key=lambda c: -c["ele"])
    return cands[0]["lon"], cands[0]["lat"], cands

PLACE_CANDS = {}   # 地名検索の候補（チェック用）

def locate_place(name, center, max_km=15):
    """地点名 → center から max_km 以内で最も近い候補 (lon, lat, title)。なければ None"""
    hits = gsi_search(name)
    exact = [h for h in hits if h[0] == name] or [h for h in hits if name in h[0]]
    near = [h for h in exact if hav(center, (h[1], h[2])) <= max_km*1000]
    PLACE_CANDS[name] = [(t, round(lo, 6), round(la, 6), round(hav(center, (lo, la))/1000, 1)) for t, lo, la in sorted(near, key=lambda h: hav(center, (h[1], h[2])))[:5]]
    if not near: return None
    t, lo, la = min(near, key=lambda h: hav(center, (h[1], h[2])))
    return lo, la, t

# ---------------------------------------------------------------- OSM
def osm_fetch(w, s, e, n):
    """OSM API map.json。範囲が大きすぎる・ノード数超過のときは4分割して統合"""
    r = S.get("https://api.openstreetmap.org/api/0.6/map.json", params={"bbox": f"{w},{s},{e},{n}"}, timeout=90)
    if r.status_code == 200: return r.json()["elements"]
    if r.status_code in (400, 509) and (e-w) > 0.01:
        mx, my = (w+e)/2, (s+n)/2; els = {}
        for b in [(w, s, mx, my), (mx, s, e, my), (w, my, mx, n), (mx, my, e, n)]:
            for el in osm_fetch(*b): els[(el["type"], el["id"])] = el
        return list(els.values())
    raise RuntimeError(f"OSM API {r.status_code}: {r.text[:200]}")

class OSM:
    HW = re.compile(r"^(path|footway|steps|track|bridleway|service|unclassified|residential)$")
    CLOSED = re.compile(r"廃道|通行止|閉鎖|立入禁止")
    def __init__(self, els, avoid=None):
        self.nodes = {e["id"]: e for e in els if e["type"] == "node"}
        self.ways = {e["id"]: e for e in els if e["type"] == "way"}
        self.adj = {}
        self.closed = [w["id"] for w in self.ways.values() if self.is_closed(w.get("tags", {}))]   # チェック用
        for w in self.ways.values():
            t = w.get("tags", {})
            if not self.HW.match(t.get("highway", "")): continue
            if t.get("access") in ("no", "private") or "disused" in t or "abandoned" in t: continue
            if self.CLOSED.search(t.get("name", "")) or (avoid and re.search(avoid, t.get("name", ""))): continue
            pen = 1.0 if t["highway"] in ("path", "footway", "steps") else (1.5 if t["highway"] == "track" else 2.0)
            for a, b in zip(w["nodes"], w["nodes"][1:]):
                if a in self.nodes and b in self.nodes:
                    d = hav(self.pos(a), self.pos(b))*pen
                    self.adj.setdefault(a, []).append((b, d, w["id"])); self.adj.setdefault(b, []).append((a, d, w["id"]))
    def pos(self, i): n = self.nodes[i]; return (n["lon"], n["lat"])
    def is_closed(self, t):
        """廃道・通行止めの道か（名前、または disused / abandoned のタグ）"""
        hw = t.get("highway") or t.get("disused:highway") or t.get("abandoned:highway") or t.get("was:highway")
        if not hw: return False
        return bool(self.CLOSED.search(t.get("name", "")) or any(k == "disused" or k == "abandoned" or k.startswith(("disused:", "abandoned:")) for k in t)
                    or t.get("highway") in ("disused", "abandoned"))
    def near(self, p):
        best = min(self.adj, key=lambda i: hav(p, self.pos(i))); return best, hav(p, self.pos(best))
    def trailheads(self):
        """登山口らしい場所（駐車場・登山口の点・「登山口」と名の付いた点）。面は重心を使う"""
        out = []
        for e in list(self.nodes.values()) + list(self.ways.values()):
            t = e.get("tags", {}); nm = t.get("name", "")
            if not (t.get("amenity") == "parking" or t.get("highway") == "trailhead" or "登山口" in nm): continue
            if e["type"] == "node": p = (e["lon"], e["lat"])
            else:
                ps = [self.pos(i) for i in e["nodes"] if i in self.nodes]
                if not ps: continue
                p = (sum(x for x, _ in ps)/len(ps), sum(y for _, y in ps)/len(ps))
            kind = "登山口" if (t.get("highway") == "trailhead" or "登山口" in nm) else "駐車場"
            out.append({"name": nm, "kind": kind, "pos": p})
        return out
    def peaks(self):
        out = []
        for n in self.nodes.values():
            t = n.get("tags", {})
            if t.get("natural") in ("peak", "volcano") and t.get("name"):
                try: ele = float(re.sub(r"[^\d.]", "", t.get("ele", "0")) or 0)
                except ValueError: ele = 0
                out.append({"name": t["name"], "ele": ele, "lon": n["lon"], "lat": n["lat"]})
        return sorted(out, key=lambda x: -x["ele"])
    def named(self):
        """名前付きの地点（ノードと、建物などの面の重心）"""
        out = []
        for e in list(self.nodes.values()) + list(self.ways.values()):
            t = e.get("tags", {}); nm = t.get("name")
            if not nm or (e["type"] == "way" and "highway" in t): continue
            if e["type"] == "node": p = (e["lon"], e["lat"])
            else:
                ps = [self.pos(i) for i in e["nodes"] if i in self.nodes]
                if not ps: continue
                p = (sum(x for x, _ in ps)/len(ps), sum(y for _, y in ps)/len(ps))
            kind = t.get("tourism") or t.get("amenity") or t.get("natural") or t.get("place") or t.get("highway") or t.get("historic") or ""
            out.append({"name": nm, "kind": kind, "lon": p[0], "lat": p[1]})
        return out

# ---------------------------------------------------------------- 地理院ベクトルタイル
def _varint(b, p):
    r = s = 0
    while True:
        c = b[p]; p += 1; r |= (c & 127) << s; s += 7
        if not c & 128: return r, p

def _zz(n): return (n >> 1) ^ -(n & 1)

def mvt_decode(b):
    import struct
    def fields(b, s, e):
        p = s
        while p < e:
            t, p = _varint(b, p); f, w = t >> 3, t & 7
            if w == 0: v, p = _varint(b, p); yield f, w, v
            elif w == 1: yield f, w, b[p:p+8]; p += 8
            elif w == 2: l, p = _varint(b, p); yield f, w, (p, p+l); p += l
            elif w == 5: yield f, w, b[p:p+4]; p += 4
            else: raise ValueError("wire type")
    def packed(rng):
        s, e = rng; out = []; p = s
        while p < e: v, p = _varint(b, p); out.append(v)
        return out
    def value(rng):
        for f, w, v in fields(b, *rng):
            if f == 1: return b[v[0]:v[1]].decode("utf-8")
            if f == 2: return struct.unpack("<f", v)[0]
            if f == 3: return struct.unpack("<d", v)[0]
            if f in (4, 5): return v
            if f == 6: return _zz(v)
            if f == 7: return bool(v)
    layers = {}
    for f, w, v in fields(b, 0, len(b)):
        if f != 3: continue
        name, keys, vals, feats, ext = "", [], [], [], 4096
        for f2, w2, v2 in fields(b, *v):
            if f2 == 1: name = b[v2[0]:v2[1]].decode()
            elif f2 == 3: keys.append(b[v2[0]:v2[1]].decode())
            elif f2 == 4: vals.append(value(v2))
            elif f2 == 5: ext = v2
            elif f2 == 2: feats.append(v2)
        out = []
        for rng in feats:
            tags, typ, geom = [], 0, []
            for f3, w3, v3 in fields(b, *rng):
                if f3 == 2: tags = packed(v3)
                elif f3 == 3: typ = v3
                elif f3 == 4: geom = packed(v3)
            props = {keys[tags[i]]: vals[tags[i+1]] for i in range(0, len(tags), 2)}
            lines, x, y, i, cur = [], 0, 0, 0, None
            while i < len(geom):
                cmd = geom[i]; i += 1; cid, cnt = cmd & 7, cmd >> 3
                for _ in range(cnt):
                    if cid == 7:
                        if cur: cur.append(cur[0])
                        continue
                    x += _zz(geom[i]); y += _zz(geom[i+1]); i += 2
                    if cid == 1: cur = [(x, y)]; lines.append(cur)
                    else: cur.append((x, y))
            out.append({"type": typ, "props": props, "lines": lines})
        layers[name] = {"extent": ext, "features": out}
    return layers

def gsi_routes(bbox, z=16):
    w, s, e, n = bbox
    x0, y1 = tile_xy(w, s, z); x1, y0 = tile_xy(e, n, z)
    jobs = [(x, y) for x in range(int(x0), int(x1)+1) for y in range(int(y0), int(y1)+1)]
    def get(xy):
        x, y = xy; r = S.get(f"{GSI}/experimental_bvmap/{z}/{x}/{y}.pbf", timeout=60)
        return (x, y, r.content if r.status_code == 200 else None)
    segs = []
    with ThreadPoolExecutor(8) as ex:
        for x, y, c in ex.map(get, jobs):
            if not c: continue
            L = mvt_decode(c).get("road")
            if not L: continue
            for f in L["features"]:
                code = f["props"].get("ftCode")
                if f["type"] != 2 or code not in (2721, 2701, 2711): continue
                for ln in f["lines"]:
                    segs.append((code, [tile_ll(x + px/L["extent"], y + py/L["extent"], z) for px, py in ln]))
    grid, NL, adj = {}, [], {}
    def node(p):
        a, b = round(p[0]*1e5), round(p[1]*1e5)
        for i in range(-2, 3):
            for j in range(-2, 3):
                k = (a+i, b+j)
                if k in grid: return grid[k]
        grid[(a, b)] = len(NL); NL.append(p); return len(NL)-1
    for code, pts in segs:
        prev = None
        for p in pts:
            nd = node(p)
            if prev is not None and prev != nd:
                d = hav(NL[prev], NL[nd]) * (1 if code == 2721 else 1.5)
                adj.setdefault(prev, []).append((nd, d, code)); adj.setdefault(nd, []).append((prev, d, code))
            prev = nd
    return NL, adj, len(jobs)

def gsi_labels(bbox, z=14):
    """地理院地図の注記（山名・地形・神社・小屋・地名）。OSMに名前のない山で地点を補う"""
    KEEP = {312: "peak", 332: "landform", 661: "worship", 673: "hut", 220: "place"}
    w, s, e, n = bbox
    x0, y1 = tile_xy(w, s, z); x1, y0 = tile_xy(e, n, z); out = []
    for x in range(int(x0), int(x1)+1):
        for y in range(int(y0), int(y1)+1):
            r = S.get(f"{GSI}/experimental_bvmap/{z}/{x}/{y}.pbf", timeout=60)
            if r.status_code != 200: continue
            L = mvt_decode(r.content).get("label")
            if not L: continue
            for f in L["features"]:
                p = f["props"]; k = KEEP.get(p.get("annoCtg"))
                if f["type"] != 1 or not k or not p.get("knj") or not f["lines"]: continue
                px, py = f["lines"][0][0]
                lo, la = tile_ll(x + px/L["extent"], y + py/L["extent"], z)
                out.append({"name": p["knj"], "kind": k, "lon": lo, "lat": la, "src": "gsi"})
    return out

# ---------------------------------------------------------------- 標高
_dem = {}
def _dem_tile(layer, z, x, y):
    k = (layer, z, x, y)
    if k not in _dem:
        r = S.get(f"{GSI}/{layer}/{z}/{x}/{y}.png", timeout=60)
        _dem[k] = Image.open(io.BytesIO(r.content)).convert("RGB").load() if r.status_code == 200 else None
    return _dem[k]

def elev(lon, lat):
    for layer, z in (("dem5a_png", 15), ("dem5b_png", 15), ("dem_png", 14)):
        fx, fy = tile_xy(lon, lat, z); tx, ty = int(fx), int(fy)
        px = _dem_tile(layer, z, tx, ty)
        if px is None: continue
        R, G, B = px[int((fx-tx)*256), int((fy-ty)*256)]
        v = R*65536 + G*256 + B
        if v == 8388608: continue
        return round((v if v < 8388608 else v-16777216)*0.01, 1)
    return None

def add_elev(pts):
    keys = set()
    for lo, la in pts:
        for layer, z in (("dem5a_png", 15), ("dem_png", 14)):
            fx, fy = tile_xy(lo, la, z); keys.add((layer, z, int(fx), int(fy)))
    with ThreadPoolExecutor(8) as ex: list(ex.map(lambda k: _dem_tile(*k), keys))  # 標高タイルを並列で先読み
    zs = [elev(lo, la) for lo, la in pts]
    nfix = sum(z is None for z in zs)
    good = [z for z in zs if z is not None]
    zs = [z if z is not None else (good[0] if good else 0) for z in zs]
    for i in range(len(zs)):  # 外れ値・欠損の補正（前後7点の中央値から80m以上ずれたら置換）
        med = statistics.median(zs[max(0, i-3):i+4])
        if abs(zs[i]-med) > 80: zs[i] = med; nfix += 1
    add_elev.last = (nfix, len(zs))
    return [(lo, la, z) for (lo, la), z in zip(pts, zs)]

def updown(p):
    return (sum(max(0, p[i][2]-p[i-1][2]) for i in range(1, len(p))),
            sum(max(0, p[i-1][2]-p[i][2]) for i in range(1, len(p))))

# ---------------------------------------------------------------- 出力
def write_kml(path, title, O, G, wps):
    def line(name, p, style):
        L = plen(p); u, d = updown(p)
        c = "\n".join(f"{x:.7f},{y:.7f},{z:.1f}" for x, y, z in p)
        return (f"<Placemark><name>{name}</name><description>距離 約{L/1000:.2f} km / 登り +{u:.0f} m・下り -{d:.0f} m</description>"
                f"<styleUrl>#{style}</styleUrl><LineString><tessellate>1</tessellate><coordinates>\n{c}\n</coordinates></LineString></Placemark>")
    pts = "\n".join(
        f"<Placemark><name>{w['name']}</name><description>スタートから約{w['km']:.2f} km / 標高 約{w['ele']:.0f} m{('／'+w['note']) if w.get('note') else ''}</description>"
        f"<styleUrl>#{'key' if i in (0, len(wps)-1) else 'wpt'}</styleUrl><Point><coordinates>{w['lon']:.6f},{w['lat']:.6f},{w['ele']:.1f}</coordinates></Point></Placemark>"
        for i, w in enumerate(wps))
    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>{title}</name>
<description>登山道：© OpenStreetMap contributors (ODbL) / 国土地理院 ベクトルタイル。標高：国土地理院 標高タイル。現地の状況・通行規制は必ず最新情報を確認してください。</description>
<Style id="osm"><LineStyle><color>ff0000ff</color><width>4</width></LineStyle></Style>
<Style id="gsi"><LineStyle><color>ffff7f00</color><width>3</width></LineStyle></Style>
<Style id="key"><IconStyle><scale>1.2</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/red-stars.png</href></Icon></IconStyle></Style>
<Style id="wpt"><IconStyle><scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href></Icon></IconStyle></Style>
<Folder><name>ルート</name>
{line('ルート（OSMの登山道）', O, 'osm')}
{line('比較：ルート（地理院地図の登山道）', G, 'gsi') if G else ''}
</Folder><Folder><name>ポイント</name>
{pts}
</Folder></Document></kml>"""
    open(path, "w", encoding="utf-8").write(kml)
    import xml.etree.ElementTree as ET; ET.parse(path)

def write_profile(path, title, O, G, wps):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fams = ["Noto Sans CJK JP", "DejaVu Sans"]
    if FONT and FONT != FONT_DEFAULT:   # 既定の Linux 以外では、見つけたフォントを登録して先頭に置く
        from matplotlib import font_manager
        try:
            font_manager.fontManager.addfont(FONT)
            fams = [font_manager.FontProperties(fname=FONT).get_name(), "DejaVu Sans"]
        except Exception:
            pass
    plt.rcParams["font.family"] = fams
    cO = [c/1000 for c in cum(O)]; zO = [p[2] for p in O]
    fig, ax = plt.subplots(figsize=(11, 4.6), dpi=130)
    ax.fill_between(cO, zO, min(zO)-50, color="#f4a6a6", alpha=.5)
    ax.plot(cO, zO, color="#d62728", lw=2, label="OSMの登山道")
    if G: ax.plot([c/1000 for c in cum(G)], [p[2] for p in G], "--", color="#1f77b4", lw=1.2, label="地理院地図の登山道")
    span = max(cO[-1], 0.1); prev, lvl = -9, 0
    for w in wps:  # 近い地点はラベルを段違いにして重ならないようにする
        lvl = lvl + 1 if (w["km"] - prev) < span*0.06 else 0; prev = w["km"]
        ax.plot(w["km"], w["ele"], "ko", ms=3)
        ax.annotate(w["name"], (w["km"], w["ele"]), fontsize=8, rotation=30, xytext=(2, 6 + 16*lvl), textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color="0.6", lw=.6) if lvl else None)
    u, d = updown(O)
    ax.set_title(f"{title}　距離 {cO[-1]:.2f} km／登り +{u:.0f} m・下り -{d:.0f} m")
    ax.set_xlabel("距離 (km)"); ax.set_ylabel("標高 (m)"); ax.grid(alpha=.3); ax.legend(loc="lower right")
    ax.set_ylim(min(zO)-50, max(zO) + (max(zO)-min(zO))*0.35 + 80)
    ax.set_xlim(-span*0.02, span*1.12)
    fig.tight_layout(); fig.savefig(path)
    fig.savefig(path.replace(".png", "_preview.png"), dpi=60); plt.close(fig)  # 確認用の縮小版

def write_map(path, title, O, G, wps, maxw=2200, maxh=1500):
    pts = [p[:2] for p in O] + [p[:2] for p in (G or [])]
    w, e = min(p[0] for p in pts), max(p[0] for p in pts); s, n = min(p[1] for p in pts), max(p[1] for p in pts)
    pw, ph = (e-w)*0.08+0.002, (n-s)*0.08+0.0015; w, e, s, n = w-pw, e+pw, s-ph, n+ph
    k = math.cos(math.radians((s+n)/2))   # 縦横比を 0.9〜2.2 に収める（細長い地図にしない）
    wm, hm = (e-w)*k, (n-s)
    if wm < 0.9*hm: d = (0.9*hm - wm)/k/2; w, e = w-d, e+d
    elif wm > 2.2*hm: d = (wm/2.2 - hm)/2; s, n = s-d, n+d
    z = 17
    while z > 10:
        x0, y0 = tile_xy(w, n, z); x1, y1 = tile_xy(e, s, z)
        if (x1-x0)*256 <= maxw and (y1-y0)*256 <= maxh: break
        z -= 1
    x0, y0 = tile_xy(w, n, z); x1, y1 = tile_xy(e, s, z)
    W, H = int((x1-x0)*256), int((y1-y0)*256)
    img = Image.new("RGB", (W, H), "white")
    def get(xy):
        tx, ty = xy; r = S.get(f"{GSI}/std/{z}/{tx}/{ty}.png", timeout=60)
        return tx, ty, (Image.open(io.BytesIO(r.content)).convert("RGB") if r.status_code == 200 else None)
    jobs = [(tx, ty) for tx in range(int(x0), int(x1)+1) for ty in range(int(y0), int(y1)+1)]
    with ThreadPoolExecutor(8) as ex:
        for tx, ty, im in ex.map(get, jobs):
            if im: img.paste(im, (int((tx-x0)*256), int((ty-y0)*256)))
    img = img.convert("RGBA"); ov = Image.new("RGBA", img.size, (0, 0, 0, 0)); dr = ImageDraw.Draw(ov)
    P = lambda p: ((tile_xy(p[0], p[1], z)[0]-x0)*256, (tile_xy(p[0], p[1], z)[1]-y0)*256)
    dr.line([P(p) for p in O], fill=(220, 0, 0, 215), width=6, joint="curve")
    if G:
        pg = [P(p) for p in G]
        for i in range(0, len(pg)-1, 2): dr.line([pg[i], pg[i+1]], fill=(0, 90, 255, 235), width=3)
    img = Image.alpha_composite(img, ov); dr = ImageDraw.Draw(img)
    f13, f15 = font(15), font(18)
    for wp in wps[1:-1]:
        x, y = P((wp["lon"], wp["lat"])); dr.ellipse([x-4, y-4, x+4, y+4], fill="black")
        tw = dr.textlength(wp["name"], font=f13); tx = x+7 if x+7+tw < W-6 else x-7-tw   # 右端からはみ出すなら左に置く
        dr.text((tx, y-20), wp["name"], font=f13, fill="black", stroke_width=3, stroke_fill="white")
    for wp, col, lab in ((wps[0], (0, 160, 0), "S"), (wps[-1], (210, 0, 0), "G")):
        x, y = P((wp["lon"], wp["lat"])); dr.ellipse([x-11, y-11, x+11, y+11], fill=col, outline="white", width=2)
        dr.text((x-5, y-12), lab, font=f15, fill="white")
        tw = dr.textlength(wp["name"], font=f15); tx = x+14 if x+14+tw < W-6 else x-14-tw
        dr.text((tx, y-12), wp["name"], font=f15, fill=col, stroke_width=3, stroke_fill="white")
    dr.rectangle([0, 0, 520, 78], fill=(255, 255, 255, 230))
    dr.text((10, 4), title, font=f15, fill="black")
    dr.line([(12, 42), (52, 42)], fill=(220, 0, 0), width=6); dr.text((60, 31), "OSMの登山道（メイン）", font=f13, fill="black")
    dr.line([(12, 64), (22, 64)], fill=(0, 90, 255), width=3); dr.line([(30, 64), (40, 64)], fill=(0, 90, 255), width=3)
    dr.text((60, 53), "地理院地図の登山道（比較）", font=f13, fill="black")
    # 縮尺バー
    mpp = 156543.03*math.cos(math.radians((s+n)/2))/2**z
    for sc in (2000, 1000, 500, 200, 100):
        if sc/mpp < W*0.25: break
    bx = sc/mpp; dr.rectangle([W-bx-30, H-48, W-20, H-8], fill=(255, 255, 255, 220))
    dr.line([(W-bx-22, H-18), (W-22, H-18)], fill="black", width=4); dr.text((W-bx-22, H-44), f"{sc} m", font=f13, fill="black")
    dr.text((10, H-24), "地図：国土地理院（地理院タイル）／登山道：© OpenStreetMap contributors, 国土地理院", font=font(12), fill=(40, 40, 40), stroke_width=2, stroke_fill="white")
    img = img.convert("RGB"); img.save(path, optimize=True)
    pv = img.copy(); pv.thumbnail((760, 760)); pv.save(path.replace(".png", "_preview.png"), optimize=True)  # 確認用の縮小版
    return z, W, H

# ---------------------------------------------------------------- 自動チェック
TRAIL_HW = ("path", "footway", "steps")
ROAD_HW = ("track", "service", "unclassified", "residential")
PLACE_SHIFT_M = 150   # 地名検索とOSMの同名地点がこれ以上離れていたらOSMを使う
BAD_PLACE_KIND = ("camp_site", "village", "hamlet", "locality", "town", "city", "suburb", "neighbourhood", "isolated_dwelling", "farm")

def build_cmd(a, **ov):
    """今の引数の一部を差し替えた、そのまま実行できるコマンド"""
    d = {"mountain": a.mountain, "start": a.start, "goal": a.goal, "start_name": a.start_name, "goal_name": a.goal_name,
         "goal_peak": a.goal_peak, "avoid": a.avoid, "title": a.title}
    d.update(ov)
    parts = ["python3", "trailkml.py"]
    for k, v in d.items():
        if v: parts += ["--" + k.replace("_", "-"), str(v)]
    if a.radius != 6.0: parts += ["--radius", str(a.radius)]
    if a.no_gsi: parts.append("--no-gsi")
    parts += ["--out", a.out]
    return " ".join(shlex.quote(x) for x in parts)

def latlon(p): return f"{p[1]:.6f},{p[0]:.6f}"

def run_checks(a, osm, path, ways, O, O2, G2, wps, picked, mt, sd, gd, efix, rep):
    """結果がおかしい可能性を調べる。[{level: 要確認|参考, code, msg, suggest}]（ルートや出力は変えない）"""
    C = []
    def add(level, code, msg, suggest=None): C.append({"level": level, "code": code, "msg": msg, "suggest": suggest})
    seglen = [hav(O2[i], O2[i+1]) for i in range(len(O2)-1)]
    total = sum(seglen) or 1
    hw = [osm.ways[w].get("tags", {}).get("highway", "") for w in ways]

    # ① ROAD_APPROACH：スタートから最初の登山道までの車道・林道
    road, k = 0.0, 0
    while k < len(hw) and hw[k] in ROAD_HW: road += seglen[k]; k += 1
    if road >= 500 or (road > 0 and road >= total*0.2):
        j = O[min(k, len(O)-1)]
        # 切り替わる地点から400m以内で、ルートの車道・林道沿い（60m以内）にある駐車場・登山口があれば、そちらを候補にする
        rd = O2[:k+1] if k > 0 else O2[:2]
        ths = [t for t in osm.trailheads() if hav(t["pos"], j[:2]) <= 400 and dist_to_line(t["pos"], rd) <= 60]
        msg = (f"スタートから最初の登山道まで車道・林道を約{road/1000:.2f}km（全体の{road/total*100:.0f}%）通っています。登山口が麓すぎる可能性があります。"
               f"車道・林道から登山道に切り替わる地点: {latlon(j)} 標高約{j[2]:.0f}m")
        if ths:
            th = min(ths, key=lambda t: (t["kind"] != "登山口", hav(t["pos"], j[:2])))
            msg += f"。その近くの{th['kind']}{('「'+th['name']+'」') if th['name'] else ''}（登山口の候補）: {latlon(th['pos'])} 標高約{elev(*th['pos']) or 0:.0f}m"
            sug = build_cmd(a, start=latlon(th["pos"]), start_name="登山口（推定）")
        else:
            msg += "（ここを登山口の候補とします）"
            sug = build_cmd(a, start=latlon(j), start_name="登山口（推定）")
        add("要確認", "ROAD_APPROACH", msg, sug)

    # ② SUMMIT_MISMATCH：終点が本当の山頂でない可能性（--goal 山頂 のとき）
    pk = picked.get("peaks_all") or []
    if "goal" in picked and picked["goal"]["by"] == "peak":
        g = picked["goal"]; why = []
        if gd > 30: why.append(f"終点が登山道から{gd:.0f}m離れています")
        top_ele = max((x["ele"] for x in pk), default=0)
        if top_ele and abs(O[-1][2] - top_ele) > 30: why.append(f"終点の標高（約{O[-1][2]:.0f}m）が山頂候補の最高標高（{top_ele:.0f}m）と{abs(O[-1][2]-top_ele):.0f}m違います")
        higher = [x for x in pk if x["ele"] > (g.get("ele") or 0) + 5 and hav((x["lon"], x["lat"]), g["pos"]) < 1500]
        if higher: why.append(f"近く（1.5km以内）に、使った山頂「{g['name']}」より高い山頂候補「{higher[0]['name']}（{higher[0]['ele']:.0f}m）」があります")
        if why:
            cands = [(x["name"], x["ele"], latlon((x["lon"], x["lat"]))) for x in pk[:5]]
            if a.goal_peak: sug = build_cmd(a, goal_peak=None)
            else:
                alt = next((x for x in pk if x["name"] != g["name"]), None)
                sug = build_cmd(a, goal_peak=alt["name"]) if alt else None
            if a.goal_peak and a.goal_peak in g["name"]:   # 山頂を名前で指定し、そのとおりに決まったときは参考にとどめる
                add("参考", "SUMMIT_MISMATCH", f"--goal-peak の指定どおり「{g['name']}」を終点にしました。" + "。".join(why) + f"。山頂候補: {cands}")
            else:
                add("要確認", "SUMMIT_MISMATCH", "終点が本当の山頂でない可能性があります。" + "。".join(why) + f"。山頂候補: {cands}", sug)

    # ③ PLACE_DOUBT：地点名の取り違え
    for label, jp in (("start", "開始点"), ("goal", "終了点")):
        p = picked.get(label)
        if not p or p["by"] == "peak": continue
        why = []
        if mt and hav(p["pos"], mt) >= 8000: why.append(f"山の位置から{hav(p['pos'], mt)/1000:.1f}km離れています")
        if p.get("kind") in BAD_PLACE_KIND: why.append(f"OSMの種類が「{p['kind']}」で、登山口らしくありません")
        if why:
            add("要確認", "PLACE_DOUBT", f"{jp}「{p['name']}」（{latlon(p['pos'])}）は取り違えの可能性があります。" + "。".join(why) + f"。他の候補: {p.get('cands') or 'なし'}")

    # ③' PLACE_SHIFT：地名検索の位置をOSMの位置に置き換えた
    for label, jp in (("start", "開始点"), ("goal", "終了点")):
        p = picked.get(label)
        if p and p.get("shift"):
            add("参考", "PLACE_SHIFT", f"{jp}「{p['name']}」は、地名検索の位置（{latlon(p['shift']['gsi'])}）がOSMの位置（{latlon(p['pos'])}）と{p['shift']['m']}m離れていたため、OSMの位置を使いました")

    # ④ DETOUR：遠回り
    straight = hav(O2[0], O2[-1])
    if straight > 0 and total/straight > 3:
        add("参考", "DETOUR", f"ルートの距離（{total/1000:.2f}km）が開始点と終了点の直線距離（{straight/1000:.2f}km）の{total/straight:.1f}倍です。遠回りしている可能性があります")

    # ⑤ GRADIENT：勾配が不自然
    if O[-1][2] > O[0][2]:
        u, _ = updown(O); gr = u/total
        if gr < 0.05 or gr > 0.6:
            add("参考", "GRADIENT", f"平均勾配（登り累計÷距離）が{gr*100:.0f}%です（目安 5〜60%）。ルートか標高を確認してください")

    # ⑥ CLOSED_NEARBY：廃道・通行止めの道が近い
    bw, be = min(p[0] for p in O2)-0.003, max(p[0] for p in O2)+0.003
    bs, bn = min(p[1] for p in O2)-0.003, max(p[1] for p in O2)+0.003
    far = [p for p in (G2 or []) if dist_to_line(p, O2) > 50]
    route_ids = set(ways)
    for wid in osm.closed:
        if wid in route_ids: continue
        w = osm.ways[wid]; t = w.get("tags", {})
        pts = [osm.pos(i) for i in w["nodes"] if i in osm.nodes]
        if len(pts) < 2 or not any(bw <= x <= be and bs <= y <= bn for x, y in pts): continue
        samp = []   # 約10mごとに点を打つ
        for i in range(1, len(pts)):
            n = max(1, int(hav(pts[i-1], pts[i])/10))
            samp += [(pts[i-1][0]+(pts[i][0]-pts[i-1][0])*j/n, pts[i-1][1]+(pts[i][1]-pts[i-1][1])*j/n) for j in range(n)]
        samp.append(pts[-1])
        nm = t.get("name") or "無名"; tagd = ",".join(k for k in t if "disused" in k or "abandoned" in k)
        label = f"「{nm}」" + (f"（{tagd}）" if tagd else "")
        near = [q for q in samp if dist_to_line(q, O2) <= 30]
        if len(near)*10 >= 100:
            mid = near[len(near)//2]
            add("要確認", "CLOSED_NEARBY", f"廃道・通行止めの道{label}が、OSMのルートから30m以内を約{len(near)*10}m並走しています（{latlon(mid)} 付近）。ルートがその道を通っていないか確認してください")
        hit = [f for f in far if min(hav(f, q) for q in samp) <= 50]
        if hit:
            mid = hit[len(hit)//2]
            add("参考", "CLOSED_NEARBY", f"地理院の線がOSMから50m以上ずれる区間（{latlon(mid)} 付近）が、廃道・通行止めの道{label}の近くを通っています。地理院の線が廃道を通っている可能性があります（メインのOSMのルートには影響しません）")

    # ⑦ ELEV_FIX：標高の補正が多い
    nfix, n = efix
    if n and nfix/n > 0.05:
        add("参考", "ELEV_FIX", f"標高の外れ値・欠損として補正した点が{nfix}/{n}点（{nfix/n*100:.0f}%）あります。登り・下りの数字は不確かです")

    # ⑧ NO_WAYPOINTS：途中の地点がない
    if len(wps) <= 2:
        add("参考", "NO_WAYPOINTS", "途中の地点情報がありません。公開されている山行記録の合目などで補うとよいです")

    # 参考：開始点・終了点が登山道から離れている（②の対象は除く）
    if sd > 30: add("参考", "SNAP_FAR", f"開始点から登山道まで{sd:.0f}m離れています")
    if gd > 30 and not any(c["code"] == "SUMMIT_MISMATCH" for c in C):
        add("参考", "SNAP_FAR", f"終了点から登山道まで{gd:.0f}m離れています（山頂の点の位置や、山頂直下の道がデータにない可能性）")
    return C

# ---------------------------------------------------------------- メイン
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mountain", default="", help="山名（地名検索で位置を決め、山頂候補の絞り込みに使う）")
    ap.add_argument("--start", required=True); ap.add_argument("--goal", required=True)
    ap.add_argument("--start-name"); ap.add_argument("--goal-name")
    ap.add_argument("--goal-peak", help="山頂候補が複数あるとき、使う山頂の名前（部分一致）")
    ap.add_argument("--avoid", help="避ける道の名前（正規表現）")
    ap.add_argument("--title"); ap.add_argument("--out", required=True)
    ap.add_argument("--radius", type=float, default=6.0, help="地点名をOSMから探すときの半径km")
    ap.add_argument("--no-gsi", action="store_true")
    a = ap.parse_args()
    rep = {"warnings": []}
    if not FONT: print("警告: 日本語フォントが見つかりません。経路図の文字が正しく表示されない可能性があります", file=sys.stderr)

    SUMMIT = ("山頂", "summit", "頂上")
    sp = parse_point(a.start); gp = None if a.goal in SUMMIT else parse_point(a.goal)
    mt = None
    if a.mountain:
        m = locate_mountain(a.mountain)
        if m:
            mt = (m[0], m[1]); rep["mountain"] = {"lon": round(m[0], 6), "lat": round(m[1], 6), "candidates": m[2][:5]}
            if len(m[2]) > 1 and m[2][1]["ele"] > m[2][0]["ele"]*0.7:
                rep["warnings"].append(f"同名の山が複数あります（標高の最も高いものを使用）: {[(c['title'], round(c['lon'], 3), round(c['lat'], 3), c['ele']) for c in m[2][:4]]}")
        else: rep["warnings"].append(f"地名検索で「{a.mountain}」が見つかりませんでした")
    center = mt or sp or gp
    if center is None:
        sys.exit("座標を指定するか、地名検索で見つかる山名（--mountain）を指定してください")
    pre = {}   # 地名検索で見つかった地点
    for label, given, p in (("start", a.start, sp), ("goal", a.goal, gp)):
        if p or given in SUMMIT: continue
        r = locate_place(given, center)
        if r:
            pre[label] = ((r[0], r[1]), r[2]); rep.setdefault("geocoded", {})[label] = [r[2], round(r[0], 6), round(r[1], 6)]
    pts = [p for p in (sp, gp) if p] + [v[0] for v in pre.values()]
    if a.goal in SUMMIT and mt: pts.append(mt)
    if len(pts) >= 2:   # 2点が決まっていれば、その範囲＋余白だけ取る
        w = min(p[0] for p in pts)-0.012; e = max(p[0] for p in pts)+0.012; s = min(p[1] for p in pts)-0.01; n = max(p[1] for p in pts)+0.01
    else:               # 地点名をOSMの中から探す必要があるときは、中心から半径 --radius km を取る
        c = pts[0] if pts else center; rad = a.radius/111
        w, e = c[0]-rad/math.cos(math.radians(c[1])), c[0]+rad/math.cos(math.radians(c[1])); s, n = c[1]-rad, c[1]+rad
    t0 = time.time(); osm = OSM(osm_fetch(w, s, e, n), a.avoid)
    rep["osm_fetch_sec"] = round(time.time()-t0, 1)
    named = osm.named()
    picked = {}   # チェック用：地点名・山頂から決めた地点

    def resolve(label, given, cand_label):
        p = parse_point(given)
        if p: return p, cand_label or ("スタート" if label == "start" else "ゴール")
        if label in pre:
            # 地名検索の施設の点は、実物から数百mずれていることがある（例：浄土平ビジターセンターは約580m西）。
            # OSMに同じ名前の地点（建物の面・点）があり、PLACE_SHIFT_M 以上離れていれば、OSMの位置を使う
            gp_, gt = pre[label]
            same = [x for x in named if x["name"] in (given, gt)]
            if same:
                o = min(same, key=lambda x: hav((x["lon"], x["lat"]), gp_)); op = (o["lon"], o["lat"]); dd = hav(op, gp_)
                if dd >= PLACE_SHIFT_M:
                    picked[label] = {"by": "osm", "name": o["name"], "kind": o["kind"], "pos": op, "cands": PLACE_CANDS.get(given, []),
                                     "shift": {"gsi": gp_, "m": round(dd)}}
                    rep["warnings"].append(f"「{given}」の地名検索の位置（{latlon(gp_)}）がOSMの「{o['name']}」（{latlon(op)}）と{dd:.0f}m離れていたため、OSMの位置を使いました")
                    return op, cand_label or o["name"]
            picked[label] = {"by": "gsi", "name": pre[label][1], "kind": "", "pos": pre[label][0], "cands": PLACE_CANDS.get(given, [])}
            return pre[label][0], cand_label or pre[label][1]
        if given in SUMMIT:
            pk = osm.peaks()
            if mt: pk = [x for x in pk if hav((x["lon"], x["lat"]), mt) < 2500] or pk   # 山名の位置の近くの山頂に絞る（隣の高い山を選ばない）
            picked["peaks_all"] = list(pk)
            if a.goal_peak: pk = [x for x in pk if a.goal_peak in x["name"]] or pk
            if not pk: sys.exit("範囲内に山頂（natural=peak）が見つかりません。終了点を座標で指定してください")
            top = pk[0]
            # 最高点の近く（1km以内）に「山頂」と名の付いた点があればそちらを優先する。
            # 例：鳥海山は「鳥海山(2236)」の点が新山と七高山の間にあり、本当の最高点は「鳥海山山頂(新山)」（標高タグなし）
            is_top = lambda nm: "山頂" in nm and not re.search(r"小屋|頂上|神社|荘|ヒュッテ|避難|トイレ", nm)   # 「月山頂上小屋」などは除く
            near_top = [x for x in pk if is_top(x["name"]) and hav((x["lon"], x["lat"]), (top["lon"], top["lat"])) < 1000]
            rep["peak_candidates"] = pk[:6] + [x for x in near_top if x not in pk[:6]]
            if near_top:
                rep["warnings"].append(f"最高点の点「{top['name']}」の近くに「{near_top[0]['name']}」があり、こちらを山頂として使いました")
                top = dict(near_top[0], ele=near_top[0].get("ele") or top["ele"])
            picked[label] = {"by": "peak", "name": top["name"], "ele": top["ele"], "pos": (top["lon"], top["lat"])}
            return (top["lon"], top["lat"]), cand_label or f"{top['name']}（{top['ele']:.0f}m）"
        hits = [x for x in named if x["name"] == given] or [x for x in named if given in x["name"]]
        if not hits: sys.exit(f"「{given}」が地名検索にもOSMの地点名にも見つかりません。座標で指定してください")
        PRI = ("trailhead", "parking", "alpine_hut", "wilderness_hut", "shelter", "peak", "information", "toilets")
        hits.sort(key=lambda h: (h["kind"] not in PRI, hav((h["lon"], h["lat"]), center)))   # 登山口・駐車場・小屋などを優先し、次に山に近い順
        rep.setdefault("name_candidates", {})[label] = [(h["name"], h["kind"], round(h["lon"], 6), round(h["lat"], 6)) for h in hits[:8]]
        if len({(round(h['lon'], 3), round(h['lat'], 3)) for h in hits}) > 1:
            rep["warnings"].append(f"「{given}」に該当する地点が複数あり、登山口・小屋などを優先し山に近い候補を使いました: {rep['name_candidates'][label][:4]}")
        picked[label] = {"by": "osm", "name": hits[0]["name"], "kind": hits[0]["kind"], "pos": (hits[0]["lon"], hits[0]["lat"]),
                         "cands": rep["name_candidates"][label][1:6]}
        return (hits[0]["lon"], hits[0]["lat"]), cand_label or hits[0]["name"]

    SP, sname = resolve("start", a.start, a.start_name)
    GP, gname = resolve("goal", a.goal, a.goal_name)
    s_id, sd, g_id, gd = snap_pair(osm.adj, osm.pos, SP, GP)
    rep.update(start=SP, goal=GP, snap_start_m=round(sd), snap_goal_m=round(gd))
    if sd > 30: rep["warnings"].append(f"開始点から登山道まで{sd:.0f}m離れています")
    if gd > 30: rep["warnings"].append(f"終了点から登山道まで{gd:.0f}m離れています（山頂の点の位置や、山頂直下の道がデータにない可能性）")
    res = dijkstra(osm.adj, s_id, g_id)
    if not res: sys.exit("OSMの登山道で開始点と終了点がつながりませんでした（--radius を広げる、地点を見直す）")
    path, ways = res
    O2 = [osm.pos(i) for i in path]
    seq = []
    for wid in ways:
        t = osm.ways[wid].get("tags", {}); nm = (t.get("name") or "無名") + (f"[{t['sac_scale']}]" if t.get("sac_scale") else "")
        if not seq or seq[-1] != nm: seq.append(nm)
    rep["osm_sequence"] = seq
    hard = [s for s in seq if re.search(r"demanding_alpine|difficult_alpine", s)]
    if hard: rep["warnings"].append(f"難易度の高い区間（sac_scale）を含みます: {hard}")

    G2 = None
    if not a.no_gsi:
        lo = [p[0] for p in O2]; la = [p[1] for p in O2]
        NL, adj, nt = gsi_routes((min(lo)-0.008, min(la)-0.006, max(lo)+0.008, max(la)+0.006))
        if adj:
            ns, _, ng, _ = snap_pair(adj, lambda i: NL[i], O2[0], O2[-1])
            r2 = dijkstra(adj, ns, ng)
            if r2:
                G2 = [NL[i] for i in r2[0]]
                dd = sorted(dist_to_line(p, O2) for p in G2)
                rep["gsi"] = {"tiles": nt, "len_m": round(plen(G2)), "snap_m": [round(hav(O2[0], NL[ns])), round(hav(O2[-1], NL[ng]))],
                              "dev_median_m": round(dd[len(dd)//2], 1), "dev_p90_m": round(dd[int(len(dd)*.9)], 1), "dev_max_m": round(dd[-1], 1)}
                # 50m以上ずれる区間の位置
                far = [p for p in G2 if dist_to_line(p, O2) > 50]
                if far: rep["gsi"]["far_points_sample"] = [(round(p[0], 5), round(p[1], 5)) for p in far[::max(1, len(far)//5)]]
            else: rep["warnings"].append("地理院の登山道では開始点と終了点がつながりませんでした")

    O = add_elev(O2); efix = add_elev.last; G = add_elev(G2) if G2 else None
    cO = cum(O)
    wps = [{"name": sname, "lon": O[0][0], "lat": O[0][1], "ele": O[0][2], "km": 0.0}]
    lo_, la_ = [p[0] for p in O2], [p[1] for p in O2]
    labels = gsi_labels((min(lo_), min(la_), max(lo_), max(la_)))
    osm_names = {x["name"] for x in named}
    for x in named + [l for l in labels if l["name"] not in osm_names]:
        if osm.CLOSED.search(x["name"]): continue          # 廃道などの名前は地点にしない
        if x["name"] in (sname, gname): continue
        p = (x["lon"], x["lat"])
        if x.get("src") == "gsi": thr = 200 if x["kind"] != "place" else 120   # 注記は文字の置き場所なので少し離れている
        else: thr = 250 if re.search(r"hut|shelter|peak|volcano|worship|place", x["kind"]) else 40
        if dist_to_line(p, O2) > thr: continue
        i = min(range(len(O)), key=lambda k: hav(p, O[k]))
        km = cO[i]/1000
        if km < 0.05 or km > cO[-1]/1000-0.05: continue
        same = [w for w in wps if hav((w["lon"], w["lat"]), p) < 40]
        if same:  # ほぼ同じ場所の地点は1つにまとめる（例：御浜小屋・御浜神社）
            if x["name"] not in same[0]["name"] and same[0]["name"].count("・") < 1: same[0]["name"] += "・" + x["name"]
            continue
        if any(w["name"] == x["name"] for w in wps): continue
        wps.append({"name": x["name"], "lon": p[0], "lat": p[1], "ele": O[i][2], "km": km, "kind": x["kind"],
                    "note": "地理院地図の注記の位置" if x.get("src") == "gsi" else ""})
    wps.append({"name": gname, "lon": O[-1][0], "lat": O[-1][1], "ele": O[-1][2], "km": cO[-1]/1000})
    wps.sort(key=lambda w: w["km"])

    title = a.title or f"{a.mountain} {sname}→{gname}".strip()
    write_kml(f"{a.out}.kml", title, O, G, wps)
    write_profile(f"{a.out}_profile.png", title, O, G, wps)
    zinfo = write_map(f"{a.out}_map.png", title, O, G, wps)
    u, d = updown(O)
    rep.update(title=title, osm_len_m=round(cO[-1]), up_m=round(u), down_m=round(d), ele_start=O[0][2], ele_goal=O[-1][2],
               ele_max=max(p[2] for p in O), waypoints=[(w["name"], round(w["km"], 2), round(w["ele"]), "gsi" if w.get("note") else "osm") for w in wps],
               map_zoom=zinfo[0], map_size=zinfo[1:], files=[f"{a.out}.kml", f"{a.out}_map.png", f"{a.out}_profile.png"])
    if G: rep["gsi"]["up_m"], rep["gsi"]["down_m"] = [round(v) for v in updown(G)]
    rep["checks"] = run_checks(a, osm, path, ways, O, O2, G2, wps, picked, mt, sd, gd, efix, rep)
    json.dump(rep, open(f"{a.out}_report.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    # 画面には要点だけ出す（詳細は _report.json）
    g = rep.get("gsi", {})
    print(f"[{title}] OSM {rep['osm_len_m']/1000:.2f}km +{rep['up_m']}m -{rep['down_m']}m 標高{rep['ele_start']:.0f}→{rep['ele_goal']:.0f}m 吸着{rep['snap_start_m']}/{rep['snap_goal_m']}m")
    if g: print(f"地理院 {g['len_m']/1000:.2f}km ずれ中央{g['dev_median_m']}m/最大{g['dev_max_m']}m" + (f" 50m超:{g['far_points_sample'][:3]}" if g.get("far_points_sample") else ""))
    if rep.get("geocoded"): print("地名検索:", rep["geocoded"])
    if rep.get("name_candidates"): print("OSM地点候補:", {k: v[:3] for k, v in rep["name_candidates"].items()})
    print("経路:", " > ".join(rep["osm_sequence"]))
    print("地点:", " / ".join(f"{n} {k}km {e}m{'*' if src == 'gsi' else ''}" for n, k, e, src in rep["waypoints"]))
    for wmsg in rep["warnings"]: print("警告:", wmsg[:300])
    for c in sorted(rep["checks"], key=lambda c: c["level"] != "要確認"):
        print(f"{c['level']}: [{c['code']}] {c['msg'][:400]}")
        if c.get("suggest"): print(f"提案: {c['suggest']}")

if __name__ == "__main__":
    main()
