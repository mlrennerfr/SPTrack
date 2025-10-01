function newtrc=shifttrc(handles,ktev,kteh,trc)
% function newtrc=shifttrc(handles,ktev,kteh,trc)
% shift correction for trc
%
% Marianne Renner - 07/10 for movtrack.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newtrc=trc;

if isempty(trc)==0
        if ktev>0 % down 
            if kteh>0 %right
                newtrc(:,3)=newtrc(:,3)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,3)=newtrc(:,3)-kteh;
            elseif kteh==0
            end
            newtrc(:,4)=newtrc(:,4)+ktev;
        elseif ktev<0 % up 
            ktev=-ktev;
            if kteh>0 %right
                newtrc(:,3)=newtrc(:,3)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,3)=newtrc(:,3)-kteh;
            elseif kteh==0
            end
             newtrc(:,4)=newtrc(:,4)-ktev;
        elseif ktev==0
            if kteh>0 %right
                newtrc(:,3)=newtrc(:,3)+kteh;
            elseif kteh<0 % left
                kteh=-kteh;
                newtrc(:,3)=newtrc(:,3)-kteh;
            elseif kteh==0
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
