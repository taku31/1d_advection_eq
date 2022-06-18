% 初期化
clearvars;

% パラメーター
Lx = 4;% 計算領域
dLx = 0.01;% 要素サイズ
dt = 0.005;
nx = Lx / dLx + 1;% 要素数
c = 1;% 移流速度
t_max = 3;
iter_max = t_max / dt;
dt_pos = 0.02;% 画像出力間隔

% クーラン数
nu = c * dt / dLx;

% 座標の生成
x = 0 : dLx : Lx;

% 配列確保
q_exact = zeros(nx, 1);
q = zeros(nx, 1);

% 初期条件の設定
for i =  1 : nx

    if abs(x(i) - 1) <= 10^-6
        q_exact(i) = 1;
        q(i) = 1;
    elseif x(i)  < 1
        q_exact(i) = 1;
        q(i) = 1;
    else
        q_exact(i) = 0;
        q(i) = 0;
    end

end

% 移流スキームの選択
method = '1st order UDS';% 1次精度風上差分法
%method = '2nd order CDS';% 2次精度中心差分法
%method = '2nd order UDS';% 2次精度風上差分法
%method = 'Lax-Friedrich';% Lax-Friedrich法
%method = 'Lax-Wendroff';% Lax-Wendroff法


% 時間発展
for iter = 1 : iter_max

    time = iter * dt;
    q_old = q;

    % 解析解
    for i = 1 : nx
        if abs(x(i) - (1 + c * time)) < 10^-6
            q_exact(i) = 1;
        elseif x(i) - (1 + c * time) < 0
            q_exact(i) = 1;
        else
            q_exact(i) = 0;
        end
    end

    % 数値解
    switch method

        case '1st order UDS'

            for i = 2 : nx

                q(i) = q_old(i) - nu * (q_old(i) - q_old(i - 1));% 1次精度風上差分

            end

        case '2nd order CDS'

            for i = 2 : nx - 1
                q(i) = q_old(i) - 0.5 * nu * (q_old(i + 1) - q_old(i - 1));% 2次精度中心差分
            end
            q(nx) = q(nx - 1);

        case '2nd order UDS'

            for i = 3 : nx
                q(i) = q_old(i) - 0.5 * nu * (3 * q_old(i) - 4 * q_old(i - 1) + q_old(i - 2));% 2次精度中心差分
            end

        case '3rd order UDS'

            for i = 3 : nx - 1
                q(i) = q_old(i) -  nu * (2 * q_old(i + 1) + 3 * q_old(i) - 6 * q_old(i - 1) + q_old(i - 2)) / 6; % 3次精度中心差分
            end

        case 'Lax-Friedrich'

            for i = 2 : nx - 1
                q(i) = 0.5 * (q_old(i + 1) + q_old(i - 1)) - 0.5 * nu * (q_old(i + 1) - q_old(i - 1));% Lax-Friedric
            end
            q(nx) = q(nx - 1);

        case 'Lax-Wendroff'

            for i = 2 : nx - 1
                q(i) = q_old(i) - 0.5 * nu * (q_old(i + 1) - q_old(i - 1)) + 0.5 * nu^2 * (q_old(i + 1) - 2 * q_old(i) + q_old(i - 1));% 2次精度中心差分
            end
            q(nx) = q(nx - 1);

        otherwise

            error("選択された手法：" + method + "はサポートしていないです。" );
    end

    % コマンドウィンドウへの出力
    txt = ['iter = ', num2str(iter), ' / ', num2str(iter_max)];
    disp(txt);

    % リアルタイム可視化
    filename = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu),'.gif'];
    fignum = 1;
    plot(x, q_exact, x, q)
    title(['time = ', num2str(time, '%.3f')]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    axis equal; axis tight; axis on;
    fig=gcf;
    fig.Color='white';
    xlim([0 4]);
    ylim([-1 2]);
    xlabel('x')
    ylabel('q')

    legtxt = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu)];
    legend('exact', legtxt,'Location','southwest')
    legend('boxoff')

    frame = getframe(fignum);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if time == dt_pos
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.05, 'Loopcount', inf);
    elseif rem(time, dt_pos) == 0
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.05, 'WriteMode', 'append');
    end

end
