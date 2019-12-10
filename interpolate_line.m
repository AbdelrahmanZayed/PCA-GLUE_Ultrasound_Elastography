function y = interpolate_line(line,y_length)

clear index

   index(1)=1;
   i=2;
   for m=1:1:y_length-1
       if line(m)~=line(m+1);
       index(i)=m;
       index(i+1)=m+1;
       i=i+2;
       end
   end
   index(end+1)=y_length;

   for u=1:2:length(index)

       if u==1
           if length(index)==2
               line(index(u):index(u+1))=linspace(line(index(u)),line(index(u+1)),index(u+1)-index(u)+1);
           else
               line(index(u):index(u+1))=linspace(line(index(u)),0.5*(line(index(u+1))+line(index(u+2))),index(u+1)-index(u)+1);
           end
       else if u~=length(index)-1
               line(index(u):index(u+1))=linspace(line(index(u-1)),0.5*(line(index(u+1))+line(index(u+2))),index(u+1)-index(u)+1);
       else
               line(index(u):index(u+1))=linspace(line(index(u-1)),line(index(u+1)),index(u+1)-index(u)+1);
           end
       end
   end

y=line;
