"4_2/post"には修士論文4.2節で使用した
重合パッチ法のポスト処理プログラムについてまとめてあります
解析プログラムで出力した変位データを主に用いて重ね合わせた変位結果やひずみ応力を算出する
また、NURBS_viewerを用いてコンター図を出力することもできる

・NURBS_input_for_s-IGA.c
    重合パッチ法ポスト処理プログラムのソースコード
    重合パッチ法解析用のNURBS_input.c（©2019年度卒業谷端さん作成NURBS_viewer）です
    .o : オブジェクトファイル

・NURBS_calc.c　（©2019年度卒業谷端さん作成NURBS_viewer）
    NURBS計算プログラムのソースコード
    .h : ヘッダーファイル、関数などについて説明してあります
    .o : オブジェクトファイル

・constant.h　（©2019年度卒業谷端さん作成NURBS_viewer）
    配列のサイズや定数の定義はここで行っています
    このファイルの変更後は"$ make clean"して"$ make"して下さい

・Makefile
    "$ make"することで実行ファイルNURBS_input_s.xが作られます
    コードの中身を書き換えたときは"$ make clean"して"$ make"して下さい

・NURBS_input_s.x
    実行ファイル
    この実行ファイルはNURBS_viewerのNURBS_inputディレクトリ内に入れておくと
    NURBS_viewerと一緒に使えて便利かと思います


・inputディレクトリ
    patch_test_integ3,6（一様引張応力問題解析）/inf_hole（無限遠中の円孔解析）/inf_crack（無限遠中の中央き裂解析）にそれぞれに分けてある

    input_for_NURBS.txt : 全てのパッチについて1つに結合したinputデータ
    input_local.txt : 全てのローカルパッチについて1つに結合したinputデータ
    Displacement.dat : 重合パッチ法解析によって取得した各コントロールポイントの変位データ
    上の3つは全て重合パッチ法解析を行うことで解析のoutputディレクトリ内に作成されます

・outputディレクトリ	
    patch_test_integ3,6（一様引張応力問題解析）/inf_hole（無限遠中の円孔解析）/inf_crack（無限遠中の中央き裂解析）にそれぞれに分けてある

    view.dat : 全ての領域におけるそれぞれのパッチの結果出力ファイル
    overlay_view.dat : ローカルパッチ領域上の重ね合わせた結果の出力ファイル
    Displacement_loc.dat : ローカルパッチ上のコントロールポイント変位
                        　 重ね合わせ結果の可視化に必要なのでここで作成する
    log_inp.dat : printfをリダイレクトしたdatファイル

    以下の8つはデータをまとめやすくするために出力したExcel等でグラフを作成するためのテキストファイルである
    disp_graph.txt : 全てのパッチ領域の変位 
    stress_vm_graph.txt : ヴォンミーゼス応力
    stress_y_graph.txt : 応力y
    over_disp_graph.txt : ローカルパッチ上の重ね合わせ変位
    over_stress_r_graph.txt : 重ね合わせ主応力
    over_stress_vm_graph.txt : 重ね合わせヴォンミーゼス応力
    over_stress_x_graph.txt : 重ね合わせ応力x
    over_stress_y_graph.txt : 重ね合わせ応力y


【コンパイル】
$ make
    Makefileを保存した場所で"$ make"する
    実行ファイルNURBS_input_s.xができる


【実行】
$ ./path/to/NURBS_input_s.x input_for_NURBS.txt Displacement.dat view.dat input_local.txt overlay_view.dat 10[xi方向分割数] 10[eta方向分割数] > ****.dat
    これはoutputディレクトリ内で行う
    view.dat : 全ての領域におけるそれぞれのパッチの結果出力ファイル
    overlay_view.dat : ローカルパッチ領域上の重ね合わせた結果の出力ファイル
    リダイレクト ">" 使ってますが直接端末に出力してもいいです
    (出力しているものが多いのでリダイレクトした方が見やすいかとは思います)
    (コアダンプするときは直接端末に出力することで、どこでコアダンプしているのか原因を見つけることができるかも)


【コマンドライン引数】
argc = 8
argv[0] = "./NURBS_input_s.x" 実行ファイル
argv[1] = "input_for_NURBS.txt" 全てのパッチについて1つに結合したinputデータ
argv[2] = "Displacement.txt" コントロールポイント変位
argv[3] = "[view].txt" NURBS_viewerに読み込ませるための全ての領域におけるそれぞれのパッチの結果出力ファイル
argv[4] = "input_local.txt" 全てのローカルパッチについて1つに結合したinputデータ
argv[5] = "[overlay_view].txt" NURBS_viewerに読み込ませるためのローカルパッチ領域上の重ね合わせた結果の出力ファイル
argv[6] = [xi方向分割数] 積分要素内の分割数
argv[7] = [eta方向分割数]


【解析結果の可視化】
可視化は"NURBS_viewer"を使います
NURBS_viewerについては2019年度卒業生谷端さんの引継データやNASを見てみてください
(念のため、私が使っていたNURBS_viewerもこのフォルダに入れておきます)
ここでは重合パッチ法解析の可視化についてのみ説明します

・全ての領域におけるそれぞれのパッチの結果出力
$ ./path/to/NURBS_view.exe input_for_NURBS.txt view.dat Displacement.dat

・ローカルパッチ領域上の重ね合わせた結果の出力
$ ./path/to/NURBS_view.exe input_local.txt overlay_view.dat Displacement_loc.dat

・NURBS_viewer表示形式の変更
NURBS_viewerの表示形式の変更はNURBS_view.c内で行うことができます
変更を行うたびに"$ make clean"と"$ make"をして下さい
NURBS_view.c内26-29行目付近ではパッチ境界やメッシュ境界の表示のデフォルトを変更できます
NURBS_view.c内600-650行目付近ではコンターバーのラベル表示のデフォルトを変更できます
どちらの項目も、NURBS_viewer実行後に変更することも出来ますが、
デフォルトを自分の好みに設定しておく方が楽だと思います
ぜひいじってみてください。