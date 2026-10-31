# mountain-trail-kml

登山道（OpenStreetMap・国土地理院）をなぞって、開始点から終了点までの最短ルートの KML・経路図・標高図を作るスクリプトです。

## 使い方

```sh
python3 trailkml.py --mountain 鳥海山 --start 鉾立 --goal 山頂 --out chokai
```

- 地点は「緯度,経度」「度分秒」「地点名」「山頂」のいずれかで指定できます。
- 出力：`<out>.kml`・`<out>_map.png`（経路図）・`<out>_profile.png`（標高図）・`<out>_report.json`
- そのほかのオプションは `python3 trailkml.py --help` を参照してください。

## 自動チェック（v2.2〜）

結果がおかしい可能性があるとき、画面に「要確認」「参考」と、そのまま実行できる修正案（「提案: python3 trailkml.py …」）を出します。`_report.json` の `checks` にも保存します。

| コード | 段階 | 内容 |
|---|---|---|
| ROAD_APPROACH | 要確認 | スタートから登山道まで車道・林道が 500m 以上、または全体の 20% 以上。近くの駐車場・登山口を候補として再実行を提案 |
| SUMMIT_MISMATCH | 要確認 | 「山頂」指定の終点が本当の山頂でない可能性（吸着 30m 超・標高差 30m 超・1.5km 以内に 5m 以上高い山頂候補） |
| PLACE_DOUBT | 要確認 | 地点名から決めた点が山から 8km 以上離れている、またはキャンプ場・集落など登山口らしくない |
| CLOSED_NEARBY | 要確認 / 参考 | 廃道・通行止めの道が OSM のルートと 30m 以内を 100m 以上並走（要確認）。地理院の線のずれが廃道付近を通る（参考） |
| DETOUR | 参考 | ルートの距離が直線距離の 3 倍を超える |
| GRADIENT | 参考 | 登りのルートで平均勾配が 5% 未満または 60% 超 |
| ELEV_FIX | 参考 | 標高の外れ値・欠損の補正が 5% 超 |
| NO_WAYPOINTS | 参考 | 途中の地点が 0 件 |
| SNAP_FAR | 参考 | 開始点・終了点が登山道から 30m 以上離れている |

## 必要なもの

- Python 3
- requests, Pillow, matplotlib
- 日本語フォント（経路図の文字に使用）。既定は `/usr/share/fonts/opentype/noto/NotoSansCJK-Medium.ttc`。無ければヒラギノなどを自動で探し、見つからなければ警告を出して続行します

## データの出典

- 登山道：© OpenStreetMap contributors (ODbL)
- 登山道・地図・標高・地名検索：国土地理院（地理院タイル・ベクトルタイル・標高タイル・地名検索）

現地の状況や通行規制は、必ず最新の情報を確認してください。
