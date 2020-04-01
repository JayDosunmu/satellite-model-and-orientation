function echo_out( out_string )

%echo_out( out_string )

global cost_params 

fid =  fopen( cost_params.out_fname,'a');

fprintf(fid,out_string);
 
fprintf(out_string);
fclose(fid);

