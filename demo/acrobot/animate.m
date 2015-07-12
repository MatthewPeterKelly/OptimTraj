function animate(t,x,P)
%animate(t,x,P)
%
%FUNCTION:
%   Animate is used to animate a system with state x at the times in t.
%
%INPUTS:
%   t = [1xM] vector of times, Must be monotonic: t(k) < t(k+1)
%   x = [NxM] matrix of states, corresponding to times in t
%   P = animation parameter struct, with fields:
%     .plotFunc = @(t,x) = function handle to create a plot
%       	t = a scalar time
%       	x = [Nx1] state vector
%     .speed = scalar multiple of time, for playback speed 
%     .figNum = (optional) figure number for plotting. Default = 1000.
%
%OUTPUTS: 
%   Animation based on data in t and x.
%

if ~isfield(P,'figNum')
    P.figNum=1000;  %Default to figure 1000
end
figure(P.figNum)


Loop_Time = 0;    %store how long has the simulation been running
T_end = t(end);   %Ending time of one step

t = t-min(t);  %Force things to start at t=0;

tic;    %Start a timer
while Loop_Time < T_end;  %Loop while the CPU time is less than the end of the simulation's time
    
    %Interpolate to get the new point:
    xNow = interp1(t',x',Loop_Time,'pchip','extrap')';
    
    %Call the plot command
   feval(P.plotFunc,Loop_Time,xNow);
   drawnow;
   
       %The next few lines pull the time step that is closest to the real time
    RealTime = toc;
    Loop_Time = RealTime*P.speed;   %Get the current time  (Taking slow motion into accunt if desired)

end

end