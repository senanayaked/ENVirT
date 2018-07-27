
%% ENVIRT_SEARCH() - Run this file to execute the ENVirT algorithm. 


function ENVIRT_SEARCH(spectrum_file, n_reads, avg_read_length, overlap_length, trim_length, M_min, M_max, p_res, p_max, L_res, L_min, L_max, n_windows, result_folder, persist_intermediates)

% [filename, path] = uigetfile('*.*', 'Select Contig Spectrum File');
% spectrum_file = fullfile(path, filename);
% contig_config_prompt = {"No. of reads", "Average Read Length", "Overlap Length", "Trim Length"};
% def_contig_config = {'', '','',''};
% contig_config = inputdlg(contig_config_prompt,"Contig Spectrum Configurations", [1, 35], def_contig_config);
% 
% n_reads = str2num(contig_config{1}); avg_read_length = str2num(contig_config{2}); 
% overlap_length = str2num(contig_config{3}); trim_length = str2num(contig_config{4});
% 
% search_config_prompt = {'M_min', 'M_max', 'p_res', 'p_max', 'L_res', 'L_min', 'L_max', 'n_windows'};
% def_search_config = {'1', '100000', '0.01', '5', '500', '500', '300500', '29' };
% 
% search_config = inputdlg(search_config_prompt, "Search Configurations", [1, 35], def_search_config);
% M_min = search_config{1};
% M_max = search_config{2};
% p_res = str2num(search_config{3});
% p_max = str2num(search_config{4});
% L_res = str2num(search_config{5});
% L_min = str2num(search_config{6});
% L_max = str2num(search_config{7});
% n_windows = str2num(search_config{8});
% 
% result_folder = uigetdir("Select Folder to Store Results");
% 
% persist_intermediates = questdlg("Would you like to save the intermediate calculation files?");

if (~exist('spectrum_file', 'var') || isempty(spectrum_file));
    error('Error, spectrum_file cannot be null');
end
if (~exist('n_reads', 'var') || isempty(n_reads)); 
    error('Error, number of reads(n_reads) cannot be null');
end
if (~exist('avg_read_length', 'var') || isempty(avg_read_length));
    error('Error, average read length(avg_read_length) cannot be null');
end
if (~exist('overlap_length', 'var') || isempty(overlap_length));
    error('Error, overlap length (overlap_length) cannot be null');
end
if ~exist('M_min', 'var') || isempty(M_min); 
    min_M = 1; 
else
    min_M = str2num(M_min);
end
if ~exist('M_max', 'var') || isempty(M_max);
    max_M = 100000; 
else
    max_M = str2num(M_max);
end
if ~exist('p_res','var') || isempty(p_res);
    p_res = 0.01;
end
if ~exist('p_max','var') || isempty(p_max); p_max = 5; end
if ~exist('L_res','var') || isempty(L_res); L_res = 500; end
if ~exist('L_min','var') || isempty(L_min); L_min = 500; end
if ~exist('L_max','var') || isempty(L_max); L_max = 300000 + L_min; end
if ~exist('n_windows','var') || isempty(n_windows); n_windows = 29; end
if ~exist('trim_length', 'var') || isempty(trim_length); trim_length = 50; end
if ~exist('result_folder', 'var') || isempty(result_folder); result_folder = 'Results'; end
if ~exist('persist_intermediates', 'var') || isempty(persist_intermediates); interm_file = false; else interm_file = persist_intermediates; end;

models = {'power law', 'exponential', 'logarithmic', 'lognormal'};

%Start parallel processing.
try
    parpool;
    mkdir(result_folder);
    C = readtable(spectrum_file, 'FileType','text', 'ReadVariableNames', false, 'ReadRowNames', false, 'Delimiter', ' ');
    C = table2array(C);
    csp_length = length(C);
    if(csp_length < trim_length);
        tail = zeros(1,(trim_length - csp_length));
        C = [C tail];
    end
    L_partition_width = 2*(L_max - L_min)/(n_windows+1);

    C = C(1:trim_length);
    C = C.*(1:length(C));
    cfg = [n_reads avg_read_length overlap_length trim_length];
    tic;
    for ab_model = 1:4    
        fprintf('%s\n',datestr(now));
        fprintf('%s\n\n',spectrum_file);

        ENVIRT_CORE(C,cfg,ab_model,min_M,max_M,p_res,p_max,L_res,L_min,L_max,L_partition_width,result_folder,interm_file);

        fprintf('%s\n\n',datestr(now));
    end
    runtime = toc;
    %End parallel processing
    fmodels = strcat(result_folder, '/model_results.txt');
    model_out = readtable(fmodels, 'FileType','text', 'ReadVariableNames', false, 'ReadRowNames', false, 'Delimiter', '\t');
    [~, inds] = sortrows(model_out(1:4, 5));
    final_out = table2array(model_out(inds(1),1:7));
    fn = strcat(result_folder, '/results_final.txt');
    fID = fopen(fn,'a');
    path = spectrum_file;
    Mf = final_out(1);
    Lf = final_out(2);
    pf = final_out(3);
    temp_m = models(1,final_out(4));
    af = temp_m{1};
    kldf = final_out(5);
    evenness = final_out(7);
    fprintf(fID,'location: %s\n richness: %i\n average genome length: %i\n d: %.5f\n a: %s\n residual error: %e\n run time: %.2fs\n evenness: %.4f\n',path,Mf,Lf,pf,af,kldf,runtime,evenness);
    fclose(fID);
        
    
    p = gcp;
    delete(p);
catch exception
    %catch exception to ensure that the parallel pool is properly closed.
    fault = exception
    p = gcp('nocreate');
    delete(p);
end
end
