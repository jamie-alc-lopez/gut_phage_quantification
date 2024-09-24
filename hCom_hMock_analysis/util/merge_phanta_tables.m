function [out_mat,master_rows] = merge_phanta_tables(tbl1,tbl2)

%This script takes in two Phanta taxa tables and produces a merged table
%with the common entries. This script only works when tbl1 and tbl2 contain
%one sample each

%Get unique set of rownames from each 
rows1 = tbl1.Taxon_Lineage_with_Names;
rows2 = tbl2.Taxon_Lineage_with_Names;
master_rows = unique([rows1; rows2]);

%Initialize empty output matrix
out_mat = zeros(length(master_rows),2);

%Loop through the combined rows 
for i = 1:size(out_mat,1)
    
    row_i = master_rows{i};

    %If that row is in tbl1 or tbl2, add the corresponding entry to
    %the merged table
    if sum(strcmp(rows1,row_i)) == 1
        out_mat(i,1) = tbl1{row_i,1};
    end
    if sum(strcmp(rows2,row_i)) == 1
        out_mat(i,2) = tbl2{row_i,1};
    end

end

end