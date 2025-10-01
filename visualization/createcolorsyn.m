function [colorm]=createcolorsyn(lastframe)
% color definitions
total=lastframe;

if lastframe<5
   colorm=[0 0 1; 0 1 1; 1 1 1; 1 1 0; 1 0 0];
else
  logical=1;
  while logical
        warning off 'all'
        C1=zeros(ceil(lastframe/8),1)';
        C2=C1;
        C3=zeros(ceil(lastframe/4),1)';
        C4=C3;
        paso=1/(lastframe/4);
        val=0;
        for ii=1:lastframe/8;
                val=val+paso;
                C1(ii)=0.5+val-0.000001;
                C2(ii)=1-val-0.000001;
                if C2(ii)<0; C2(ii)=0; end
                if C1(ii)<0; C1(ii)=0; end
                if C2(ii)>1; C2(ii)=1; end
                if C1(ii)>1; C1(ii)=1; end
        end
        val=0;
        for ii=1:lastframe/4;
                val=val+paso;
                C3(ii)=val-0.000001;
                C4(ii)=1-val-0.000001;
                if C3(ii)<0; C3(ii)=0; end
                if C4(ii)<0; C4(ii)=0; end
                if C3(ii)>1; C3(ii)=1; end
                if C4(ii)>1; C4(ii)=1; end
        end
        unos=ones(ceil(lastframe/4),1);
        ceros1=zeros(ceil(lastframe/8),1);
        ceros2=zeros(ceil(lastframe/4),1);
        colorm=[ceros1 ceros1 C1'; ceros2 C3' unos; C3' unos C4'; unos C4' ceros2; C2' ceros1 ceros1];
        clear C1 C2 C3 C4
        if size(colorm,1)>total-1
            logical=0;
            break
        else
            lastframe=lastframe+2;
        end
  end % while
end

clear ceros1 ceros2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
