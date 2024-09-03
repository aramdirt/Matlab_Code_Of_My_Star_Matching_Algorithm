
function [table] = bubble_sort(table,table_size,key_idx,ascending)

for i=1:1:table_size
    for j=1:1:table_size-i
        exchange=0;
        if ascending~=0
            if table(j,key_idx) > table(j+1,key_idx)          
                exchange=1;
            end
        else
            if table(j,key_idx) < table(j+1,key_idx)          
                exchange=1;
            end
        end
        
        if exchange==1
            tmp=table(j,:);
            table(j,:)=table(j+1,:);
            table(j+1,:)=tmp;
        end
    end
end
