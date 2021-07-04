# 環境設定方法とコードの説明
## 環境設定方法
---
1. Pythonのダウンロードとインストール
    [Pythonホームページ](https://www.python.org/)
    インストール時にAdd Python 3.9 to PATHにチェックを入れてください
    ![](https://www.javadrive.jp/python/install/img/p1-9.png)
    [参考ページ](https://www.javadrive.jp/python/install/index1.html)

2. VSCodeのダウンロードとインストール
    [VSCode ダウンロードページ](https://code.visualstudio.com/Download)

3. Python 拡張機能のインストール
   ![](https://docs.microsoft.com/ja-jp/learn/language/python-install-vscode/media/visual-studio-code-extensions-install.png)
   正常にインストールされると.pyの拡張子ファイルで右上に実行ボタンが表示されるはずです
   ![](https://i.imgur.com/5vixus9.png)
   出ない場合は，一度拡張機能をアンインストールしてからもう一度インストールしてください

4. Pythonのライブラリのインストール
    VSCodeを開き画面上部にあるターミナルをクリックして
    新しいターミナルを適当な場所に作成してください(コマンドプロンプトでも多分大丈夫です)
    ![](https://i.imgur.com/hXezxEB.png)
    1. pip自体のアップグレード
        コマンドラインで以下のコードを入力して実行してください
        ```
        py -m pip install -U pip
        ```
    2. NumPyのインストール
        コマンドラインで以下のコードを入力して実行してください
        ```
        py -m pip install numpy
        ```
    3. pandasのインストール
        コマンドラインで以下のコードを入力して実行してください
        ```
        py -m pip install pandas
        ```
    4. Matplotlibのインストール
        コマンドラインで以下のコードを入力して実行してください
        ```
        py -m pip install matplotlib
        ```
    以上で環境設定完了ですお疲れさまでした
    [参考ページ](https://gammasoft.jp/blog/install-numpy-pandas-matplotlib-by-pip/)

## コードの説明
---
```
b_spline
    b_spline_function.py        :   関数のモジュール，このファイルが実行ファイルと同じフォルダにないと正常に実行できません
    b_spline_main.py            :   基本的なIGAのプログラムです
    b_spline_surface.py         :   surfaceでのIGAプログラム，輪講本の例題を再現しています
```
以下```b_spline_main.py```について説明します
行数(line)はデフォルトの値です
    
- line 9~
  CP(コントロールポイント)の座標入力，座標は任意

- line 20
  次数入力

- line 24~
  アフィン変換(基本的には変更しなくてよい)

- line 44
  0連続となるコントロールポイント番号を入力
  2次以上はできないし，最大1点までしか入力できません
  (ファンクション適当に作った)
  コントロールポイント番号は0から始まることに注意してください
  エラーメッセージが出たらとりあえずこの行列を
  ```
  make_C0_CP = np.array([])
  ```
  に変更したら上手くいくかもしれません

- line 47~
  knot insertion Aはコントロールポイント番号を入力してその位置でノットインサーションを行います
  通常ノットインサーションの入力値は挿入したいノットの値ですが，コントロールポイントが多いといちいちノットの値を計算するのが面倒なので作りました
  例えばコントロールポイントの総数が8個で
  ```
  knot = [0. 0. 0. 0.2 0.4 0.6 0.8 0.8 1. 1. 1.]
  ```
  のとき，コントロールポイント番号2を入力すると
  ```
  new_knot_position = np.array([2])
  ```
  となり，
  ```
  knot = 0.3
  ```
  が挿入された場合と同じ結果になります
  コントロール番号が0と最後の番号(この例では7)ではこのファンクションは使えません
  よくわからない場合は次のknot insertion Bで直接ノットの値を入力してください

- line 51~
  knot insertion Bは通常のノットインサーションです
  以下のように挿入したいノットを入力してください
  例：
  ```
  insert_knot = np.array([1/3, 2/3])
  ```

- line 59~
  オーダーエレベーションを行います
  オーダーエレベーションのファンクション内では，
  特定の位置でknot insertion →
  ベジェ曲線に分割し，分割後のコントロールポイント番号を取得 →
  各ベジェ曲線でオーダーエレベーション →
  knot insertionしたノットをremoval knot(knot insertionの逆の操作)
  の手順で次数をあげています
  詳しくはNURBS bookを参照してください
  removal knotについては詳しく書かれていなかったのですが，
  knot insertionの
  $$
  \overline{B} = {\bf T}B
  $$
  の式から($\overline{B}$と${\bf T}$が既知でコントロールポイントを減らしたい)
  $$
  B = {\bf T}^{-1}\overline{B}
  $$
  となります，ただし逆行列を求めるためには${\bf T}$が正則行列である必要があるのですが，${\bf T}$は正則行列ではないので疑似逆行列という概念を用いて計算しました
  詳しくは```b_spline_function.py```を見てみてください(超読みにくい)

以上で説明は終了です
line 155~とかところどころax4ありとかなしとか書いてますが
これはオーダーエレベーションしたときに元の曲線と一致してるか確認したかっただけなので，特に気にしなくていいです．

気になる場合は
```
# ax4なしの時
```
と書いてる行のコメントアウトを外して
```
# ax4ありの時
```
と書いてる行を消すかコメントアウトにして
```
# oeder elevation 前のスプライン(黒) 比較用
#---------------------------------------------------------------#

#---------------------------------------------------------------#
```
の間を全部消すかコメントアウトにすれば消えると思います
