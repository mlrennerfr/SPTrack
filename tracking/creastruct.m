function [fit]=creastruct(trcfilename)
%function [fit]=creastruct(trcfilename)
% from a file in trc formats, creates the structure fit to save .traj file
% Marianne Renner mar 09 SPTrack v4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=dlmread(trcfilename);
size_M=length(M(:,1));
nrotrc=max(M(:,2));
if size_M~=0
    ligne=1;
    numero=1;
    prev_frame=M(1,2);
    num_seg=1;
    fit.spot(numero).nb_points=1;
    fit.spot(numero).nb_segments=1;
    fit.spot(numero).segment(num_seg).coordinates=[M(ligne,3:4) M(ligne,2)];
    fit.spot(numero).segment(num_seg).length=1;
    for ligne=2:size_M
        
        if M(ligne,1)==numero
            fit.spot(numero).nb_points=fit.spot(numero).nb_points+1;
            if M(ligne,2)==prev_frame+1
                fit.spot(numero).segment(num_seg).coordinates=[fit.spot(numero).segment(num_seg).coordinates; [M(ligne,3:4) M(ligne,2)]];
                fit.spot(numero).segment(num_seg).length=fit.spot(numero).segment(num_seg).length+1;
                prev_frame=M(ligne,2);
            else
                fit.spot(numero).nb_segments=fit.spot(numero).nb_segments+1;
                num_seg=num_seg+1;
                fit.spot(numero).segment(num_seg).coordinates=[M(ligne,3:4) M(ligne,2)];
                fit.spot(numero).segment(num_seg).length=1;
                prev_frame=M(ligne,2);
            end
        else
            numero=numero+1;
            num_seg=1;
            fit.spot(numero).nb_points=1;
            fit.spot(numero).nb_segments=1;
            fit.spot(numero).segment(num_seg).coordinates=[M(ligne,3:4) M(ligne,2)];
            fit.spot(numero).segment(num_seg).length=1;
            prev_frame=M(ligne,2);
        end
    end
    fit.nb_spots=M(size_M,1);
    fit.methode='';
end


% end of file