
%all loops start from 2 because 1 is for null complex or environment
% write in text file for Bertini
fileID = fopen(filnam ,'w');
fprintf(fileID, 'CONFIG \n \nEND; \n \nINPUT \n \nvariable_group ');
q = sym('q',[num_spec,1]);

for i = 1:num_spec-1
    fprintf(fileID, string(q(i))+', ');
end
i = num_spec;
fprintf(fileID, string(q(i))+';');

%write function list
fprintf(fileID,'\nfunction ')
r = sym('r',[num_spec,1]);
for i = 1:num_spec-1
    fprintf(fileID, string(r(i))+', ');
end
i = num_spec;
fprintf(fileID, string(r(i))+';');
fprintf(fileID,'\n\n');

for i = 1:num_spec
    fprintf(fileID,string(r(i))+' = '+string(MAK(i))+';\n');
end
fprintf(fileID,'\nEND;')