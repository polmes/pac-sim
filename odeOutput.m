function status = odeOutput(t,y,flag)
    global oT oY;

    if nargin < 3 || isempty(flag)
        oT = [oT t];
        oY = [oY y];
    else
        switch flag
            case 'init'
%                 fprintf('start\n');
%                 oT = t(1);
%                 oY = y;
            case 'done'
%                 fprintf('done\n');
        end
    end
    
    status = 0;
end
