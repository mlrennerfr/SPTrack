function newimage=shiftimage(handles,ktev,kteh,image,Ydim,Xdim)
% function newimage=shiftimage(handles,ktev,kteh,image,Ydim,Xdim)
% shift images
% Marianne Renner 09/10 SPTrack v_4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newimage=zeros(Ydim,Xdim); %att background
if isempty(image)==0
        if ktev>0 % down 
            if kteh>0 %right
                newimage(1+ktev:Ydim,1+kteh:Xdim)=image(1:Ydim-ktev,1:Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(1+ktev:Ydim,1:Xdim-kteh)=image(1:Ydim-ktev,1+kteh:Xdim);
            elseif kteh==0
                newimage(1+ktev:Ydim,:)=image(1:Ydim-ktev,:);
            end
        elseif ktev<0 % up 
            ktev=-ktev;
            if kteh>0 %right
                newimage(1:Ydim-ktev,1+kteh:Xdim)=image(1+ktev:Ydim,1:Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(1:Ydim-ktev,1:Xdim-kteh)=image(1+ktev:Ydim,1+kteh:Xdim);
            elseif kteh==0
                newimage(1:Ydim-ktev,:)=image(1+ktev:Ydim,:);
            end
        elseif ktev==0
            if kteh>0 %right
                newimage(:,1+kteh:Xdim)=image(:,1:Xdim-kteh);
            elseif kteh<0 % left
                kteh=-kteh;
                newimage(:,1:Xdim-kteh)=image(:,1+kteh:Xdim);
            elseif kteh==0
                newimage=image;
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
