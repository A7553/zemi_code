tic;
%% 
%canvas.mを実行して文字やイラストをかく．
%名前を"draw.png"で保存する．
Io = imread("zemi.png");
imwrite(Io,'image3d.bmp','bmp');

%画像を拡大する
I = imresize(Io, 3.0, 'nearest');
I = rgb2gray(I);

%二値化
BW = imbinarize(I);
%エッジ検出
BW_c = edge(BW,'canny');
%表示
hold on
imshow(BW_c)
hold off
title('Canny Filter');
axis equal

%% 
%画像サイズの取得
[x,y] = size(BW_c);
%Xは三角形分割の点の配列．
X = [];
%Cは連結リストの配列．Xのidxで連結を表現している．ex)C([1,2])・・・X[1]とX[2]の座標をつなぐ
%Cが三角分割を行う際の境界線になる．
C = [];

%探索したかどうかを保存する配列（訪れたら１にする）
seen = zeros(x,y);
%移動８方向のベクトル
step = containers.Map({0,1,2,3,4,5,6,7},{[-1,-1],[0,-1],[1,-1],[1,0],[1,1],[0,1],[-1,1],[-1,0]});

%点を表示する
figure(2);
%現在の点の総数．Cの連結リストを作る時に必要．
cnt = 0;

%全探索
for j=1:y
    for i=1:x
        if(BW_c(i,j)==true && seen(i,j)==0) %訪れていない輪郭の場合，探索を開始する
            
            hold on
            %スタート点はxで表示
            plot(i, j, 'x');
            seen(i,j) = 1;
            %スタート点は後で使うので保存する
            start = [i j];
            start_cnt = cnt + 1;
            
            prev = 1; %前の進行方向(prev->nowの向き)
            next = 1; %次の進行方向(now->nextの向き)
            now = start;
            
            %一時的なメモリ.最終的にはX,Cに結合する
            x_tmp = [];
            c_tmp = [];
            %スタート点を格納する
            x_tmp = cat(1,x_tmp,start);

            while(next~=-1) %次の点がなくなるまで探索を続ける
                %調べるのは，前回の進行方向から半時計回りに90度から
                seach = mod(prev+5,8);
                %step変数の初期化
                di = 0;
                dj = 0;
                %見つかったかどうかの判定
                find = false;
                %探索は８方向のうち，４方向のみで十分
                for k = 0:4
                    if(find == false) %見つかっていない場合
                        %反時計回りにヒトマス進める
                        seach = mod(seach+1,8);
                        d = step(seach);
                        %xの進むstep
                        di = d(1,1);
                        %yの進むstep
                        dj = d(1,2);
                        if(0<now(1,1)+di && now(1,1)+di<=x && 0<now(1,2)+dj && now(1,2)+dj<=y)
                            if(BW_c(now(1,1)+di, now(1,2)+dj)==1 && seen(now(1,1)+di, now(1,2)+dj)==0) %輪郭の探索が画像の中に収まっているかどうかの判定
                            find = true;
                            next = seach;
                            end
                            seen(now(1,1)+di, now(1,2)+dj) = 1;
                        end
                       

                    else %次の点が見つかった場合(72行目の分岐)
                        seach = mod(seach+1,8);
                        %輪郭が二重になっていた場合に内側の点を消す処理(92〜106)
                        if(mod(next,2)==1) %次の点が奇数の場合
                            if(mod(next+2,8)==seach)
                                d = step(seach);
                                di = d(1,1);
                                dj = d(1,2);
                                seen(now(1,1)+di, now(1,2)+dj) = 1;
                            end
                        else %次の点が偶数の場合
                            if(mod(next+1,8)==seach)
                                d = step(seach);
                                di = d(1,1);
                                dj = d(1,2);
                                seen(now(1,1)+di, now(1,2)+dj) = 1;
                            end
                        end
                        
                    end 

                end %for文(71行目)終了

                if(find)

                    d = step(next);
                    di = d(1,1);
                    dj = d(1,2);
                    now = [now(1,1)+di, now(1,2)+dj];
                    %直線上の点列(next=prev)は重要ではないのでX,Cに保存しない．一方で，曲線上の点列の場合(next~=prev)は保存する
                    %データの削減を行っている
                    if(next ~= prev)
                        cnt = cnt + 1;
                        %点の表示
                        plot(now(1,1), now(1,2), 'o');
                        x_tmp = cat(1,x_tmp,now);
                        c_tmp = cat(1,c_tmp,[cnt cnt+1]);
                    end
                    
                    prev = next; 
                   
                else %次の点がない場合
                    cnt = cnt + 1;
                    %start点をつなぐ
                    c_tmp = cat(1,c_tmp,[cnt start_cnt]);
                    %while文を抜けるために-1に設定する
                    next = -1;
                end   
                
            end %while終了

            %一時的なメモリx_tmp,c_tmpをX,Cに結合する
            X = cat(1,X,x_tmp);
            C = cat(1,C,c_tmp);

        end
        
    end
end

figure(2);
            
hold on
%境界線を赤色で表示している
plot(X(C'),X(C'+size(X,1)),'-r','LineWidth',2)
hold off

figure(3);
hold off
axis equal
%% 

%三角形分割を行う
%Figure3を拡大すれば見える
dt = delaunayTriangulation(X,C);
io = isInterior(dt);
figure(4);


%patch('Faces',dt(io,:),'Vertices',dt.Points,'FaceColor','w')

tr = triangulation(dt(io,:),dt.Points); triplot(tr)

patch('Faces',tr.ConnectivityList,'Vertices',tr.Points,'FaceColor','w')

%z成分(0)を追加する
Z = zeros(size(tr.Points,1),1);
%三次元モデルの下の部分
tr2.Points = cat(2,tr.Points,Z);
tr2.ConnectivityList = tr.ConnectivityList;

%z成分(100)を追加する
Z = Z + 100.0;
%三次元モデルの上の部分
tr3.Points = cat(2,tr.Points,Z);
tr3.ConnectivityList = tr.ConnectivityList + size(tr.Points,1);

tr_all.Points = cat(1,tr2.Points,tr3.Points);
tr_all.ConnectivityList = cat(1,tr2.ConnectivityList,tr3.ConnectivityList);

disp(size(tr.Points,1));
disp(y);
axis equal


%C2,C3で側面の三角形分割を行っている．上下の点をうまく繋いでいるので点は増えていない
C2 = cat(2,C,C(:,1)+size(tr.Points,1));
C3 = cat(2,C+size(tr.Points,1),C(:,2));

%C2,C3を全体の三角形分割に追加する
tr_all.ConnectivityList = cat(1,tr_all.ConnectivityList,C2);
tr_all.ConnectivityList = cat(1,tr_all.ConnectivityList,C3);
patch('Faces',tr_all.ConnectivityList,'Vertices',tr_all.Points,'FaceColor','w')

%もともと実装されているstlwriteではなく，FileExchangeからダウンロードした関数である
%FileExchangeとは，ユーザが作ったコードを共有するサイトである
stlwrite('figure_ab.stl',tr_all.ConnectivityList,tr_all.Points);

%stlファイルを書き出す
Letter = stlread('figure_ab.stl');

%モデルや始点の設定
model.Vertices = Letter.Points;
model.Faces    = Letter.ConnectivityList;
figure(5), clf
patch(model, 'FaceColor', [0.8 0.8 1.0], ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.15)
view([20 20])
lightangle(-45,70)

axis equal
xlim([0 x])
ylim([0 y])
zlim([-200 700])

view([20 20])

toc;