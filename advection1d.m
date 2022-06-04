% ������
clearvars;

% �p�����[�^�[
Lx = 4;% �v�Z�̈�
dLx = 0.1;% �v�f�T�C�Y
dt = 0.05;
nx = Lx / dLx + 1;% �v�f��
c = 1;% �ڗ����x
t_max = 3;
iter_max = t_max / dt;
dt_pos = 0.2;% �摜�o�͊Ԋu

% �N�[������
nu = c * dt / dLx;

% ���W�̐���
x = 0 : dLx : Lx;

% �z��m��
q_exact = zeros(nx, 1);
q = zeros(nx, 1);

% ���������̐ݒ�
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

% �ڗ��X�L�[���̑I��
method = '1st order UDS';% 1�����x���㍷���@
%method = '2nd order CDS';% 2�����x���S�����@

% ���Ԕ��W
for iter = 1 : iter_max

    time = iter * dt;
    q_old = q;

    % ��͉�
    for i = 1 : nx
        if abs(x(i) - (1 + c * time)) < 10^-6
            q_exact(i) = 1;
        elseif x(i) - (1 + c * time) < 0
            q_exact(i) = 1;
        else
            q_exact(i) = 0;
        end
    end

    % ���l��
    for i = 2 : nx - 1 % ���E�v�f�͌Œ�Ƃ��Čv�Z�͔�΂��B
        switch method

            case '1st order UDS'

                q(i) = q_old(i) - nu * (q_old(i) - q_old(i - 1));% 1�����x���㍷��


            case '2nd order CDS'

                q(i) = q_old(i) - 0.5 * nu * (q_old(i + 1) - q_old(i - 1));% 2�����x���S����

            otherwise

                error("�I�����ꂽ��@�F" + method + "�̓T�|�[�g���Ă��Ȃ��ł��B" + ...
                    "1st order UDS, 2nd order CDS����I�����Ă��������B");
        end
    end

    % �R�}���h�E�B���h�E�ւ̏o��
    txt = ['iter = ', num2str(iter), ' / ', num2str(iter_max)];
    disp(txt);

    % ���A���^�C������
    filename = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu),'.gif'];
    fignum = 1;
    plot(x, q_exact, x, q)
    title(['time = ', num2str(time, '%.3f')]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    axis equal; axis tight; axis on;
    fig=gcf;
    fig.Color='white';
    xlim([0 4]);
    ylim([0 1.5]);
    xlabel('x')
    ylabel('q')

    legtxt = [method, ', dt = ', num2str(dt), ', dx = ', num2str(dLx),', c = ', num2str(c),', nu = ', num2str(nu)];
    legend('exact', legtxt,'Location','southwest')
    legend('boxoff')

    frame = getframe(fignum);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if time == dt_pos
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.2, 'Loopcount', inf);
    elseif rem(time, dt_pos) == 0
        imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.2, 'WriteMode', 'append');
    end

end