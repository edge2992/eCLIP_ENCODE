## やりたいこと

## 検証

- [x] foldseekのTMscoreの行列は対象行列ではないらしいので、どうなっているかを確認して、必要であれば対象行列にすることを考える。 -> metrics_foldseek.py
- タンパク質ごとにスコアの分布を見る
- nullが入っていた場所がどこかを確認する
- score同士の分布を可視化して特徴を掴む
    - [x] アライメントスコアとfoldseekのTMscoreを比較する -> metrics_foldseek.py
    - stringdbの複数のスコアを比較する
    - stringdbのスコアとkeywordのスコアを比較する (同じkeywordで複合体を形成しているタンパク質ペアがあった時にscoreが高くあってほしい)
    - foldseekのTMscoreとkeywordのスコアを比較する (ドメインが同じであればTMscoreが高くあってほしい)


## 考察

### 配列アライメントスコアとしてのbitscoreと構造アライメントスコアTMscore

>bitscoreは偶然の一致で見つかるシーケンスデータベースに必要なサイズ, log2でスケーリングされる
>経験即では50以上のbitスコアが良いらしい
https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html

### TODO: 
- 閾値を50あたりにして分布を比較してみる
- homologyスコアをstringdbから取得して, blastpの結果と比較する
