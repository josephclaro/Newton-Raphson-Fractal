%% Chaos Project - Fractals
% Joseph Claro - 21317623
% Dynamical Systems

%% NEWTON-RAPHSON FRACTAL

%% Setup
clf

a = 7; % 
b = 8/5; %
c = 6i; %
d = 8; % polynomial coefficients
g = 0; %-1; % (a + bx + cx2 + dx3 + gx4 + hx5 + jx6 = 0)
h = 0; %-2; %
j = 0; %2; %

% Polynomial
f = @(x) 5.8.*x.^7+9; %@(x) a + b.*x + c.*x.^2 + d.*x.^3 + g.*x.^4 + h.*x.^5 + j.*x.^6;
% Derivative
df = @(x) 7*5.8.*x.^6; %@(x) b + 2*c.*x + 3*d.*x.^2 + 4*g.*x.^3 + 5*h.*x.^4 + 6*j.*x.^5;

% Newton-Raphson Iteration
nraph = @(x) x - (f(x)./df(x));

fig = 1;

%% Identifying roots
% Sparse mesh for time-saving

T = 2; B = -2; % define grid region
L = -2; R = 2; %

testPts = makeGrid(3,-3,-3,3,10); % top, bottom, left, right, # of pts down/across

testTol = 1e-6; % tight tolerance for accuracy
rootsFound = NRmethod(nraph, testPts, testTol, 1000); % Do Newton-Raphson to find roots
rootsFound = reshape(rootsFound,[],1); % A big vector of all the 'landing sites' after doing N-R

rootsSorted = sort(rootsFound); % Sort roots by abs. value, to enable plucking out unique roots

roots = rootsSorted(1); % initialise list of unique roots as the first point landed on
rootToCheck = rootsSorted(1); % " comparison root " " " ...

for j = 1:length(rootsSorted)
    if (abs(rootsSorted(j) - rootToCheck) > testTol) % if diff between current point and already-found root > tol...
        roots = [roots;rootsSorted(j)]; % add to list of 'unique' roots (if < tol, root is considered too similar)
        rootToCheck = rootsSorted(j); % comparison root is now this value
    end   
end

figure(fig)
clf
plot(real(roots),imag(roots),"mo");
title("Newton-Raphson Iteration - Result", "Interpreter","latex");
xlim([L R]); ylim([B T]);
xlabel("Re(z)","interpreter", "latex"); ylabel("Im(z)","Interpreter","latex","Rotation",pi/2);
xline(0); yline(0); grid on;
hold on;

% Stay on same fig for next plot

%% Example on one point

z = zeros(20,1); % vector of N-R iterations
z(1) = -1+1i; % initial value

for n = 2:(length(z))
    z(n) = nraph(z(n-1)); % apply Newton-Raphson to prev. point
end
figure(fig)
hold on;
plot(real(z),imag(z),"bo-");
plot(real(z(1)),imag(z(1)),"r*","markerSize",15);
title("Newton-Raphson Iteration - Roots / Example on Single Point", "Interpreter","latex");
xlim([-2 2]); ylim([-2 2]);
xlabel("Re(z)","interpreter", "latex"); ylabel("Im(z)","Interpreter","latex","Rotation",pi/2);
xline(0); yline(0); grid on;

shg

fig=fig+1;
%% Fractal-plotting Grid stuff

N = 1500; % number of points per row/col

T = 2; B = -2; % define grid region (top/bottom/left/right)
L = -2; R = 2; %

gridPts = makeGrid(T,B,L,R,N);

%% The actual Fractal
colours = ["#2e4045"; "#83adb5"; "#c7bbc9"; "#5e3c58"; "#bfb5b2"; "#b45f06";"#f7da7c";"#7edbb3";"#d690f9"];
%colours = ["#70B596"; "#66023c"; "c"; "#5e3c58"; "#bfb5b2"; "#b45f06";"#f7da7c";"#7edbb3";"#d690f9"];

% #70B596 oxidized copper
% #66023c tyrian purple

tol = 1e-4; % largest stepsize to consider a root reached

tic % start timer
% Do Newton-Raphson method from a given starting point(s)
fractRoots = NRmethod(nraph, gridPts, tol, 500); % points, tol, maxIters

colourBins = zeros(0.5*N^2,length(colours)); %initialise list(s), as matrix, to plot as same colour

colourIndex = ones(length(colours),1); % used to add points to above colour bins

for i = 1:N
    for j = 1:N
        for c = 1:length(roots)
        % if root is... add to respective colour vector
            if (abs(fractRoots(i,j) - roots(c)) < tol)
                colourBins(colourIndex(c),c) = gridPts(i,j);
                colourIndex(c) = colourIndex(c) + 1;
                break
            end
        end
    end
end

figure(fig)
clf
hold on;
title("Newton-Raphson Fractal for Pippi hehehe", "Interpreter","latex");
xlim([L R]); ylim([B T]);
xlabel("Re(z)","interpreter", "latex"); ylabel("Im(z)","Interpreter","latex","Rotation",pi/2);
xline(0); yline(0); grid on;

for c = 1:length(roots)
    plot(real(colourBins(:,c)),imag(colourBins(:,c)),".","Color",colours(c),"MarkerSize",2);
end

plot(real(roots),imag(roots),"m.","MarkerSize",10);

shg
toc
fig=fig+1;

%% Functions

function grid = makeGrid(top,bottom,left,right,N)

    meshX = linspace(left,right,N);
    meshY = linspace(top,bottom,N);
    
    gridPts = zeros(length(meshY),length(meshX));
    
    for m = 1:length(meshX) % assign real parts (- to +)
        gridPts(:,m) = gridPts(:,m) + meshX(m)*ones(length(meshY),1);
    end
    
    for n = 1:length(meshY) % assign imaginary parts (+ to -)
        gridPts(n,:) = gridPts(n,:) + meshY(n)*1i*ones(1,length(meshY));
    end
    grid = gridPts;
end

function result = NRmethod(nraph, pts, tol, maxIters)

    maxStep=1;
    iters = 0;
    x_old = pts;
    
    while (maxStep > tol || iters < maxIters)
        x_new = nraph(x_old); % x_new = x_old - f/df
      
        maxStep = max(abs(x_new - x_old),[],"all"); % find largest step out of all values
        x_old = x_new;
    
        iters=iters+1;
    end
    result = x_old;
end
