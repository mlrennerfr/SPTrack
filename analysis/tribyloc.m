function cut=tribyloc(tri,perival)

if perival==1 % peri=syn
    count=1;
    nb=tri.nrosegm;
    if nb==1
        cut=tri;
    else
        cut.segment(count).data=[];
        for order=1:tri.nrosegm-1
            if isempty(tri.segment(order).data)==0 & isempty(tri.segment(order+1).data)==0
                if tri.segment(order).data(1,6)~=0 & tri.segment(order).data(1,6)==-tri.segment(order+1).data(1,6)
                    cut.segment(count).data=[cut.segment(count).data; tri.segment(order).data];
                else
                    cut.segment(count).data=[cut.segment(count).data; tri.segment(order).data];
                    count=count+1;
                    cut.segment(count).data=[];
                end
            end
        end
        cut.segment(count).data=[cut.segment(count).data; tri.segment(nb).data];
        cut.data=tri.data;
        cut.nrosegm=count;
    end
    
else % peri=extra
    count=1;
    if tri.nrosegm==1
        cut=tri;
    else
        cut.segment(count).data=[];
        for order=1:tri.nrosegm-1
            if isempty(tri.segment(order).data)==0 & isempty(tri.segment(order+1).data)==0
                if tri.segment(order).data(1,6)<=0 & tri.segment(order+1).data(1,6)<=0
                    cut.segment(count).data=[cut.segment(count).data; tri.segment(order).data];
                else
                    cut.segment(count).data=[cut.segment(count).data; tri.segment(order).data];
                    count=count+1;
                    cut.segment(count).data=[];
                end
            end
        end
        cut.segment(count).data=[cut.segment(count).data; tri.segment(tri.nrosegm).data];
        cut.data=tri.data;
        cut.nrosegm=count;
    end
end %perival