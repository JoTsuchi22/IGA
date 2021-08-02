"4_2/analysis"には修士論文4.2節で使用した
重合パッチ法の解析プログラムについてまとめてあります

・s_IGA_array_patch_test_integ3.c
    重合パッチ法プログラム（ガウス積分の積分点3x3）
    4_2節"重合パッチ法を用いた解析"では基本的にこのプログラムを使用している

・s_IGA_array_patch_test_integ6.c
    重合パッチ法プログラム（ガウス積分の積分点6x6）
    4.2.1.3節"積分点数の違いによる調査"で使用した

これらの解析プログラムでは、全コントロールポイントの変位データを取得することが目的である
重ね合わせ結果やひずみ応力結果についてはその後のポスト処理(4_2/postに示す)で計算・可視化する

重合パッチ法では今のところzarusoba_viewerは使わず、
NURBS_viewerを使うようなプログラム構成にしている

3次元の重合パッチ法を行う場合、
NURBS_viewerを3次元化する
or
3次元zarusoba_viewerに即した形式のファイルを出力するようにコードを修正する
必要がある

・inputディレクトリ
    patch_test_integ3,6（一様引張応力問題解析）/inf_hole（無限遠中の円孔解析）/inf_crack（無限遠中の中央き裂解析）にそれぞれに分けてある
    
    input_glo_**.txt : グローバルパッチのinputデータ
    input_loc_**.txt : ローカルパッチのinputデータ
    input_format.txt : inputデータのフォーマット
                    　 重合パッチ法で用いるinputデータはグローバル・ローカル共に通常のIGAと全く同じものである（2019年度修了長島彩華さんの引継データやNASも参考にして下さい）

・outputディレクトリ
    patch_test_integ3,6（一様引張応力問題解析）/inf_hole（無限遠中の円孔解析）/inf_crack（無限遠中の中央き裂解析）にそれぞれに分けてある
	
    checkAns：コントロールポイントと積分点上での解析結果出力用のディレクトリです。
 		        基本の確認はこのディレクトリ内のデータで大まかに解析チェックしてます。
	colored_point：zarusoba_viewerで点郡に色をつけて可視化するためのデータがあるディレクトリ。
		        （今は空だがプログラムのコメントアウトをはずせば使うことができるはず）
	Gauss_stress：zarusoba_viewerで積分点に色を付けて可視化するためのディレクトリ
  	mesh_net：メッシュ図とコントロールネット図を作るためのデータが入るディレクトリ。
		        （今は空だがプログラムのコメントアウトをはずせば使うことができるはず）
	new_zarusoba：zarusoba_viewerで変位を可視化するためのディレクトリ
		        （今は空だがプログラムのコメントアウトをはずせば使うことができるはず）
    NURBS：NURBSで補間した解析モデルと変位のデータ等が入るディレクトリです。
	shapefunc：形状関数の形をグラフ化するためのデータが入るディレクトリ。

    coord_data.txt : グローバル・ローカル内の全てのコントロールポイント番号とそれらの座標
                     各コントロールポイントに与えた座標についてExcelなどで簡単に確認することができる
    Displacement.dat : 各コントロールポイント座標の変位データ
                       重合パッチ法解析ではこのデータを用いてポスト処理(ひずみや応力の結果算出)を行うので、
                       このデータに特に注意するようにして下さい(-nanやinfなどが出ていないか、など)
    input_for_NURBS.txt : 重合パッチ法解析に使用した複数のinputデータ(グローバルとローカル)を1つのinputデータにしたもの
                        　NURBS_viewerによる可視化の際に使用する
    input_local.txt : 解析に使用しローカルパッチのinputデータ
                    　NURBS_viewerによる可視化の際に使用する
    test.dat : printfをリダイレクトしたdatファイル


【コンパイル】
$ gcc -g -o s_IGA_array_patch_test_integ*.x path/to/s_IGA_array_patch_test_integ*.c -lm -Wall -mcmodel=large


【実行前の作業】
$ cd path/to/[outputディレクトリを作りたい場所]
    解析結果をまとめたい場所に移動

$ mkdir output_*
    output_*ディレクトリに解析結果をまとめる
    出力データはたくさんあるので上のコマンドで解析結果出力ディレクトリを作った方がいいと思います

$ cd output_*
    outputディレクトリに移動

$ mkdir checkAns colored_point Gauss_stress mesh_net new_zarusoba NURBS shapefunc
    このコマンドで以下に示す7つのディレクトリを作成できます
    checkAns / colored_point / Gauss_stress / mesh_net / new_zarusoba / NURBS / shapefunc
    これを忘れると解析がコアダンプします、漏れがないよう作成してください


【解析実行】
$ ./path/to/s_IGA_array_patch_test_integ*.x path/to/input_glo_*.txt path/to/input_loc_*.txt > ****.dat
    これはoutputディレクトリ内で行う
    input_glo_*.txt : グローバルパッチのinputデータ
    input_loc_*.txt : ローカルパッチのinputデータ
    リダイレクト ">" 使ってますが直接端末に出力してもいいです
    (出力しているものが多いのでリダイレクトした方が見やすいかとは思います)
    (もしコアダンプしてしまう場合は、リダイレクトせず直接端末に出力することで、どこでコアダンプしているのか原因を見つけることができるかも、試してみてください)


【コマンドライン引数】
argc = 3~ (ローカルパッチの総数によって変わる)
argv[0] = "./[実行ファイル].x"
argv[1] = "[input_global].txt"
argv[2] = "[input_local1].txt"
これ以降はローカルパッチが複数あるとき
argv[3] = "[input_local2].txt"
			・
			・
			・