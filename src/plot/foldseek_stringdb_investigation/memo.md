# やりたいこと

## 検証

- [x] foldseekのTMscoreの行列は対象行列ではないらしいので、どうなっているかを確認して、必要であれば対象行列にすることを考える。 -> metrics_foldseek.py
- [x] タンパク質ごとにスコアの分布を見る
- [x] nullが入っていた場所がどこかを確認する
- score同士の分布を可視化して特徴を掴む
  - [x] アライメントスコアとfoldseekのTMscoreを比較する -> metrics_foldseek.py
  - [x] stringdbの複数のスコアを比較する -> metrics_stringdb.py
  - [x] stringdbのスコアとkeywordのスコアを比較する (同じkeywordで複合体を形成しているタンパク質ペアがあった時にscoreが高くあってほしい) -> keyword_stringdb.py
  - [x] foldseekのTMscoreとkeywordのスコアを比較する (ドメインが同じであればTMscoreが高くあってほしい) -> keyword_foldseek.py

## 考察

### 配列アライメントスコアとしてのbitscoreと構造アライメントスコアTMscore

>bitscoreは偶然の一致で見つかるシーケンスデータベースに必要なサイズ, log2でスケーリングされる
>経験即では50以上のbitスコアが良いらしい
[参考](https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html)

### TMScoreはクエリのタンパク質の長さが短ければ、高く出る傾向がある

統計に影響するのでmax, min, avgの処理を追加する

### TODO

- [x] 閾値を50あたりにして分布を比較してみる
- [x] homologyスコアをstringdbから取得して, blastpの結果と比較する
