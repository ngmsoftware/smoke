clear();
clc();
close('all');

figure('position',[100 100 800 800]);

sat = @(t, s)(t.*(t<s)+s.*(t>=s));

makeMovie = 1;

S = configSimlation('0','periodicCircular');

gridSizeX = S.gridSizeX;
gridSizeY = S.gridSizeY;
X = S.X;
Y = S.Y;
nu = S.nu;
k = S.k;
dt = S.dt;
N = S.N;
srcPos = [0.5 0.5];
vx = S.vx;
vy = S.vy;

p = zeros(gridSizeY, gridSizeX);


if makeMovie
    vid = VideoWriter('output/test.avi');
    open(vid);
end

cla();
hold('on');
H = imagesc(p,[0 1]);
delta = 8;
Q = quiver(gridSizeX*(X(1:delta:end,1:delta:end)+1)/2, gridSizeY*(Y(1:delta:end,1:delta:end)+1)/2, vx(1:delta:end,1:delta:end) , vy(1:delta:end,1:delta:end),'color','w');
axis([0 gridSizeX 0 gridSizeY])

set(Q,'HitTest','off');
set(get(Q,'Parent'),'HitTest','off')
set(H, 'ButtonDownFcn','t=0;');

t = 0.0;
for iter=1:N
    t = t+dt;

    
    % SOURCE and FIELD configuration
    %
    
    
    %pos = S.srcPos(t);
    %p( 1+fix(pos(1)*S.gridSizeY), 1+fix(pos(2)*S.gridSizeX) ) = S.srcAmp(t)/5;
    
%     
%     p(5, 1+fix(0.5*S.gridSizeX)) = 3;
%     vy(round(5:end/8), 1+fix(0.5*S.gridSizeX)) = 0.4;
%     vx(5:10, 1+fix(0.5*S.gridSizeX)) = 0.2*cos(2*t);
% 
%     p(end-5, 1+fix(0.5*S.gridSizeX)) = 3;
%     vy(end-round(5:end/8), 1+fix(0.5*S.gridSizeX)) = -0.4;
%     vx(end-round(5:10), 1+fix(0.5*S.gridSizeX)) = 0.2*sin(2*t);
% 
%     p(1+fix(0.5*S.gridSizeX),5) = 3;
%     vx(1+fix(0.5*S.gridSizeX),round(5:end/8)) = 0.4;
%     vy(1+fix(0.5*S.gridSizeX),5:10) = 0.2*cos(2*t);
% 
%     p(1+fix(0.5*S.gridSizeX),end-5) = 3;
%     vx(1+fix(0.5*S.gridSizeX),end-round(5:end/8)) = -0.4;
%     vy(1+fix(0.5*S.gridSizeX),end-round(5:10)) = 0.2*sin(2*t);
    
    
    
     p(end/4,end/4) = 3;
     p(end/4+(-2:2),end/4+(-2:2)) = conv2(p(end/4+(-2:2),end/4+(-2:2)),[0 1 0; 1 1 1; 0 1 0]/5,'same');
     vx(end/4,end/4) = 6*(0.3+0.3*sin(10*t));
     vy(end/4,end/4) = 6*(0.3+0.3*cos(10*t));

     vx(end/2,end/2) = 1.0;
     vy(end/2,end/2) = 1.0;
    
    
    
    
    vx = conv2(vx,[0 1 0; 1 1 1; 0 1 0]/5,'same');
    vy = conv2(vy,[0 1 0; 1 1 1; 0 1 0]/5,'same');
    
    
    [p, vx, vy] = navierStokesStep(k, p, nu, vx, vy, dt, srcPos);
    
    div = -divergence(vx,vy);
    idiv = reshape(poicalc(div(:),1,1,size(div,1),size(div,2)),size(div,1),size(div,2));
    [gx, gy] = gradient(idiv);
    vx = vx-gx;
    vy = vy-gy;    
    
    vx = sat(vx,3);
    vy = sat(vy,3);
    
    if mod(iter,12)==0
        set(H,'CData',p);
        set(Q,'udata',vx(1:delta:end,1:delta:end) ,'vdata',vy(1:delta:end,1:delta:end));
        drawnow();
        
        if makeMovie
            currFrame = getframe;
            writeVideo(vid,currFrame);    
        end
        
    end

end

if makeMovie
    close(vid);
end