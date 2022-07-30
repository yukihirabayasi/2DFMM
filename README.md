2D Fast Multipole Method(FMM) 
====

渦法への実装は"Install in Vortex Solver (ver C)"を読んでください 

※1: Fortranは大規模な構造体を動的配列で用いた場合に、サイズやランクなどの付加情報がC言語に比べて非常に大きくなり実行困難になることや、モートンオーダーとの相性がやや悪いことからC言語へ移植しFortran版は廃止とした。これに伴いFortran版の動作保証はされていない。  
※2: 見込み角によって直接計算するべきと判定されたものはd_biot.f90ではなく、fmm.cのdirect_evaluateで計算されている点に注意。  
※3: mirrorには対応してません。  

# Files
./main_solver/testFMM.f90 ... 渦法ソルバーに実装する前にs\_\*\*\*\*\_step.datのみでTree法が正しく計算できているか確認するためのテストソルバー  
./debug.py ... 渦要素位置や生成された格子、モートンオーダーの実装などを確認できる  
./multipole_test.f90 ... S2Mが正しく行われているのか確認できる  
./read_s_file.f90 ... s\_\*\*\*\*\_step.datを読み込むプログラム  
./m_tree_bindc.f90 ... C言語で書かれたTree法の関数をFortranで書かれた渦法ソルバーで呼び出せるようにする  
./morton_order.c ... モートンオーダーに関する必須のライブラリ  
./quadtree.f90 ... Tree構造に関する必須のソルバー  
./fmm.f90 ... 生成されたTree構造上で速度場を計算するための必須のソルバー  

# Functions
渦法ソルバー側
+ Make_tree / ツリー構造を作成し、S2MからM2Mまでを行う
+ biot_tree / M2Sを実行し、速度場を算出する

Tree法本体
編集中

# Test Run (ver Fortran)
~~・src内に含まれるコードで構成されたメイン機能  
　Makeでコンパイルを行い、testFMM.outを実行するのみ  
・debug.py  
　testFMM.f90を実行することで生成されるvor.datが必要となる~~  

# Test Run (ver C) 
$ icc -c fmm.c quadtree.c morton_order.c  
$ ifort fmm.o quadtree.o morton_order.o testFMM.f90

# Install in Vortex Solver (ver C)
ver-c/sample_implementは渦法にTree法を試験的に実装されたもので、これと同じようにすれば動く。また、e_velo_mpi.f90の編集(以下3〜5)はファイルごと置き換えることも有効。
1. 渦法ソルバー/src/treeにこのリポジトリのver-c/src内fmm.c fmm.h quadtree.c quadtree.h morton_order.c morton_order.h m_tree_bindc.f90を配置する
2. cal_cond/fmm_param.datを作成
3. Tree法を用いて誘起速度の計算を行うソルバー(e_velo_mpi.f90など)でtree/m_tree_bindc.f90をinclude
4. biot_bを行う以前(宣言文の直下など)にcall make_tree(nvor_b,vor_b,cir_b,cor_b)を記述
5. biot_bをbiot_treeに書き換える
6. Tree法のソルバーをコンパイル`mpicc -c tree/morton_order.c tree/fmm.c tree/quadtree.c`
7. 渦法のソルバーをコンパイル`mpif90 fmm.o quadtree.o morton_order.o a_main.f90`
8. 実行`./a.out`、`a_main.exe`など

## トラブルシューティング(詳細なプログラムの説明)
渦法ソルバー側で実行する必要があるのはm_tree_bindc.f90のmake_tree、biot_treeの2つのみなので、プログラムとしてはmake_treeにblob渦の情報を渡し、biot_treeに算出したいx,y座標を渡すだけで実行できる。ただし、make_treeとbiot_treeはC言語で書かれているのでこれをfortran側で再定義するためにm_tree_bindc(モジュール)を読み込む必要がある。またmake_tree、biot_tree本体であるプログラム(morton_order.c、quadtree.c、fmm.c)をsrc内に配置し、cal_cond内に設定ファイルであるfmm_cond.datを追加する必要がある。
