function InvPend03_Animate(t,xx,u,u_min,u_max,J,fname)

% animation credits to Renato Coral, UnB
% link to thesis: https://repositorio.unb.br/bitstream/10482/34629/1/2018_RenatoCoralSampaio.pdf

v = VideoWriter(fname,'Motion JPEG AVI');
open(v);

maxpos = 1.2*max(abs(xx(1,:)));
maxth  = 1.2*max(abs(xx(2,:)));
maxJ   = 1.2*max(abs(J));
minJ   = 1.2*min(abs(J));

% size of the cart
H = 0.0; % horizon line
L = 0.4*maxpos; % length of the pendulum for ploting

% Inverted Pendulum Problem

rSmpl = 2;

niter = length(t);
for i=1:rSmpl:niter
    
    fprintf('Frame %d out of %d \n',i,niter);
    
    figure(1); clf; 
    l = 2;
    
    p  = xx(1,i); 
    th = xx(2,i);
    
    set(gcf, 'pos',[0 40 800 1000]);
    whitebg([1.0 1.0 1.0])
    set(gcf,'Color',[1 1 1])
    subplot(6,1,[1 2]); hold all
    % line([-obj.state_lower_limits(1) obj.state_upper_limits(1)], [H H], 'color', 'k', 'Linewidth',2.5); hold on;
    % line([obj.state_lower_limits(1) obj.state_upper_limits(1)], [H H], 'color', 'k', 'Linewidth',2.5); hold on;
    
    title(['Simulation Time: ' num2str(t(i)) 's']);
    plot(p,H,'ks','MarkerSize',20,'Linewidth',3); % posicao do carrinho
    
    xB = p - L*sin(th);
    yB = H + L*cos(th);
    
    plot([p;xB],[H;yB],'r-o','Linewidth',l); % pendulo
    
    grid on;
    xlim([-maxpos maxpos]);
    ylim([-maxpos maxpos]);
    axis square
    
    subplot(6,1,3);
    plot(t(1:i),xx(1,1:i),'linewidth',1,'color','k'); hold on
    yline(0,'r--','linewidth',1.5),
    ylim([-5 maxpos])
    xlim([0 max(t)])
    grid on;
    ylabel('p','FontSize',13)
    
    subplot(6,1,4);
    plot(t(1:i),xx(2,1:i),'linewidth',1,'color','k'); hold on
    yline(0,'r--','linewidth',1.5),
    ylim([-maxth maxth])
    xlim([0 max(t)])
    grid on;
    ylabel('th','FontSize',13)
    
    subplot(6,1,5);
    plot(t(1:i),u(1:i),'k','linewidth',2), hold on
    yline(u_max,'r--','linewidth',1.5),
    yline(0,'r--','linewidth',1.5),
    yline(u_min,'r--','linewidth',1.5),
    xlim([0 max(t)])
    ylabel('u','FontSize',13);
    grid on;
        
    % cost function
    subplot(6,1,6);
    semilogy(t(1:i),J(1:i),'linewidth',1.5,'color','k','linestyle','-','marker','.');hold on
    yline(0,'r--','linewidth',1.5),
    ylim([minJ*1e-1 maxJ*1e1])
    xlim([0 max(t)])
    grid on;
    ylabel('Cost function','FontSize',13);
    
    frame = getframe(gcf);
    writeVideo(v,frame);    
end

close(v);