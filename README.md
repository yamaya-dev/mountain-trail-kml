# mountain-trail-kml

登山道（OpenStreetMap・国土地理院）をなぞって、開始点から終了点までの最短ルートの KML・経路図・標高図を作るスクリプトです。

## 使い方

```sh
python3 trailkml.py --mountain 鳥海山 --start 鉾立 --goal 山頂 --out chokai
```

- 地点は「緯度,経度」「度分秒」「地点名」「山頂」のいずれかで指定できます。
- 出力：`<out>.kml`・`<out>_map.png`（経路図）・`<out>_profile.png`（標高図）・`<out>_report.json`
- そのほかのオプションは `python3 trailkml.py --help` を参照してください。

## 必要なもの

- Python 3
- requests, Pillow, matplotlib
- Noto Sans CJK フォント（経路図の文字に使用。パスは `/usr/share/fonts/opentype/noto/NotoSansCJK-Medium.ttc`）

## データの出典

- 登山道：© OpenStreetMap contributors (ODbL)
- 登山道・地図・標高・地名検索：国土地理院（地理院タイル・ベクトルタイル・標高タイル・地名検索）

現地の状況や通行規制は、必ず最新の情報を確認してください。
