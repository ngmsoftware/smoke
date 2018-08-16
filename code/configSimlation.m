function S = configSimlation(vectorFieldName, sourceName)

S.gridSizeX = 200;
S.gridSizeY = 200;

[X, Y] = meshgrid(linspace(-1,1,S.gridSizeX), linspace(-1,1,S.gridSizeY));
S.X = X;
S.Y = Y;

S.vx = zeros(S.gridSizeY);
S.vy = zeros(S.gridSizeX);

% S.nu = 0.001;
% S.k = .0001;
% S.dt = 0.0025;
% S.N = 16000;
% S.t = 0.0;


S.nu = 0.01;
S.k = .0001;
S.dt = 0.0015;
S.N = 16000;
S.t = 0.0;


S.srcPos = [0.5 0.5];


%vx((X.^2+Y.^2)<0.1) = 0;
%vy((X.^2+Y.^2)<0.1) = 0;
 



switch vectorFieldName
    case '0'

    case 'XY'
        S.vx = -.4*sign(X).*abs(X).^.5;
        S.vy = -.4*sign(Y).*abs(Y).^.5;

        z = sqrt(X.^2+Y.^2)>0.25;
        
        S.vx = S.vx.*z;
        S.vy = S.vy.*z;
        
    case 'simple'
        
        S.vy = .5*ones(S.gridSizeX);
        S.vx = .5*ones(S.gridSizeY);
        

    case 'circle'
        for i=1:S.gridSizeY
            for j=1:S.gridSizeX
                x = X(i,j);
                y = Y(i,j);

                V1 = -[0 1; -1.0 0]*[y; x] + [0.0; 0.0];
                S.vx(i,j) = .4*V1(2);
                S.vy(i,j) = .4*V1(1);
                
            end
        end
    
        
    case 'square'
        for i=1:S.gridSizeY
            for j=1:S.gridSizeX
                x = X(i,j);
                y = Y(i,j);

                S.vx(i,j) = 0.8*y.^7;
                S.vy(i,j) = -0.8*x.^7;

                n = sqrt(S.vx(i,j).^2+S.vy(i,j).^2);
                
                S.vx(i,j) = S.vx(i,j)/n;
                S.vy(i,j) = S.vy(i,j)/n;
                
            end
        end
        
        
        
    case 'sincos'
        
        for i=1:S.gridSizeY
            for j=1:S.gridSizeX
                x = X(i,j);
                y = Y(i,j);

                 V1(1) = 0.1*cos(8*x);
                 V1(2) = -0.1*sin(8*y+pi/2);
                 S.vx(i,j) = 80*V1(2)/20;
                 S.vy(i,j) = 80*V1(1)/20;
            end
        end

        
        
        
        
    case 'fillipov'
        
        for i=1:S.gridSizeY
            for j=1:S.gridSizeX
                x = X(i,j);
                y = Y(i,j);

                A = [0.7578 -1.9796; 1.7454 -0.3350];
                b = [0.1005; -2.1600];
                a = 6.2759;
                v = [-0.1582; 1.8467];
                xx = 20*[x; y];
                V = A*xx + a*b*sign(v'*xx);
                S.vx(i,j) = 0.015*V(1);
                S.vy(i,j) = 0.015*V(2); 
            end
        end
end


switch sourceName
   
    
    case 'simple'
        S.srcPos = @(t)[0.5 0.5];
        S.srcAmp = @(t)1;
        
        
    case 'randPos'
        S.srcPos = @(t)rand(1,2).*(mod(t,0.05)<1.01*S.dt);
        S.srcAmp = @(t)80;
        
        
        
    case 'periodicFixed'
        S.srcPos = @(t)[0.5 0.5];
        S.srcAmp = @(t)80.0*sin(5*t).^48;

        
    case 'periodicCircular'
        S.srcPos = @(t)[0.4*sin(4*t)+0.5, 0.4*cos(4*t)+0.5];
        S.srcAmp = @(t)60.0*sin(20*t).^48;
        
    case 'periodicCircular2'
        S.srcPos = @(t)[0.4*sin(43*t)+0.5, 0.4*cos(43*t)+0.5];
        S.srcAmp = @(t)6.0*cos(.2*t).^48;
end
