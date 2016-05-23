function ENVIRT_CORE(config_spectrum,cfg,ab_model,min_M,max_M,p_res,p_max,L_res,L_min,L_max,L_partition_width, result_folder, interm_file)

tt = tic;

if (min_M > 100000)
    min_res = 5;
end
if (min_M <= 100000)
    min_res = 4;
end
if (min_M <= 10000)
    min_res = 3;
end
if (min_M <= 1000)
    min_res = 2;
end
if (min_M <= 100)
    min_res = 1;
end

models = ['power law', 'exponential', 'logarithmic', 'lognormal'];

min_M_int = round(min_M/10^(min_res-1));

max_res = round(min(log10(max_M)-1,5));
max_M_int = round(max_M/10^(max_res-1));

ab = [1,2,3,4]; b=ab_model;

n_partitions = (L_max - L_min)/(L_partition_width/2)-1;
remainder = rem((L_max - L_min),L_partition_width/2);

if(remainder > 0)
    n_partitions = (L_max - L_min - remainder)/(L_partition_width/2);
end

ga_results = inf(max_res*n_partitions + 2, 10);
ga_file = strcat(result_folder,'/ga_results_',num2str(ab(b)),'.txt');
csp = config_spectrum;
for i = min_res:max_res

    M_res = 10^(i-1);

    res = [M_res,L_res,p_res];
    M_LB = 1; M_UB = 100;
    if (i == min_res)
        M_LB = min_M_int;
    end
    if (i == max_res)
        M_UB = max_M_int;
    end
    p_LB = 1; p_UB = round(p_max/p_res);
    
    for j = 1:n_partitions

        L_LB = (L_min+(j-1)*L_partition_width/2)/L_res; 
        L_UB = (L_min+(j+1)*L_partition_width/2)/L_res; 
        LB = [M_LB,L_LB,p_LB]; UB = [M_UB,L_UB,p_UB];

        [x,fval,ef,runTime] = ga_iteration(i,j,ab(b),res,LB,UB,n_partitions,csp,cfg);
        
        if (interm_file);
            fn = strcat(result_folder,'/results.txt');
            fID = fopen(fn,'a');
            fprintf(fID,'%i\t%i\t%.5f\t%i\t%e\t%i\t%.4f\n',x(1),x(2),x(3),models(ab(b)),fval,ef,runTime);
            fclose(fID);
        end; 
        row_n = (i-1)*n_partitions + j;
        ga_results(row_n,:) = [LB(1)*M_res,UB(1)*M_res,LB(2)*L_res,UB(2)*L_res,x(1),x(2),x(3),ab(b),fval,ef];     

    end

    fprintf('%s\n','-------------------------------------------------------------');

end
% Final refinement stage
[min_row,~] = find(ga_results == min(ga_results(:,9)));
min_row = min(min_row);

L_LB = max(ga_results(min_row,6) - 1000,L_min);
L_UB = min(ga_results(min_row,6) + 1000,L_max);
L_res = 1;

M_LB = max(ga_results(min_row,5) - 1000,1);
M_UB = min(ga_results(min_row,5) + 1000,max_M);
M_res = 1;

p_res = p_res/10;
p_LB = max((ga_results(min_row,7)-0.25)/p_res,1);
p_UB = (ga_results(min_row,7)+0.25)/p_res;

res = [M_res,L_res,p_res];
LB = [M_LB,L_LB,p_LB]; UB = [M_UB,L_UB,p_UB];

[x,fval,ef,~] = ga_iteration(max_res+1,1,ab(b),res,LB,UB,n_partitions,csp, cfg);

row_n = max_res*n_partitions + 1;
ga_results(row_n,:) = [LB(1)*M_res,UB(1)*M_res,LB(2)*L_res,UB(2)*L_res,x(1),x(2),x(3),ab(b),fval,ef];
dlmwrite(ga_file,ga_results);

[min_row,~] = find(ga_results == min(ga_results(:,9)));
min_row = min(min_row);

Mf = ga_results(min_row,5);
Lf = ga_results(min_row,6);


boundry = ((Mf == M_LB)||(Mf == M_UB)||(Lf == L_LB)||(Lf == L_UB));
not_max_it = 1;
n_it = 0;

while (boundry && not_max_it)

    n_it = n_it + 1;

    M_LB = max(ga_results(min_row,5) - 5000,1);
    M_UB = min(ga_results(min_row,5) + 5000,max_M);

    L_LB = max(ga_results(min_row,6) - 5000,L_min);
    L_UB = min(ga_results(min_row,6) + 5000,L_max);

    p_LB = max((ga_results(min_row,7)-0.25)/p_res,1);
    p_UB = (ga_results(min_row,7)+0.25)/p_res;

    res = [M_res,L_res,p_res];
    LB = [M_LB,L_LB,p_LB]; UB = [M_UB,L_UB,p_UB];

    [x,fval,ef,~] = ga_iteration(max_res+1,1+n_it,ab(b),res,LB,UB,n_partitions,csp, cfg);

    row_n = max_res*n_partitions + 1 + n_it;
    ga_results(row_n,:) = [LB(1)*M_res,UB(1)*M_res,LB(2)*L_res,UB(2)*L_res,x(1),x(2),x(3),ab(b),fval,ef];
    dlmwrite(ga_file,ga_results);

    [min_row,~] = find(ga_results == min(ga_results(:,9)));
    min_row = min(min_row);

    Mf = ga_results(min_row,5);
    Lf = ga_results(min_row,6);

    boundry = ((Mf == M_LB)||(Mf == M_UB)||(Lf == L_LB)||(Lf == L_UB));
    not_max_it = (n_it < 3);

end


[min_row,~] = find(ga_results == min(ga_results(:,9)));
min_row = min(min_row);

Mf = ga_results(min_row,5);
Lf = ga_results(min_row,6);
pf = ga_results(min_row,7);
af = ga_results(min_row,8);
kldf = ga_results(min_row,9);
evenness = calculateEvenness(af, pf, Mf);

rt = toc(tt);
fn = strcat(result_folder, '/model_results.txt');
fID = fopen(fn,'a');
fprintf(fID,'%i\t%i\t%.5f\t%i\t%e\t%.4f\t%.4f\n',Mf,Lf,pf,af,kldf,rt,evenness);%'location: %s\n richness: %i\n average genome length: %i\n d: %.5f\n a: %s\n residual error: %e\n run time: %.2fs\n evenness: %.4f\n',path,Mf,Lf,pf,models(af),kldf,rt,evenness);
fclose(fID);
fprintf('%i\t%i\t%.5f\t%s\t%e\t%.4f\t%.4f\n',Mf,Lf,pf,models(af),kldf,rt,evenness);
dlmwrite(ga_file,ga_results);
end

%The cost function definition%

function y = cost_fn(x, cfg, csp, abm, run)

R = cfg(1); r = cfg(2); o = cfg(3); 

C = csp;

ab = abm;
M_res = run(1);L_res = run(2);p_res = run(3);

M = round(x(1))*M_res;
L = round(x(2))*L_res;
d = x(3)*p_res;

E = zeros(1,length(C));
Z = zeros(1,length(C));

H = zeros(M,1);
d = round(d*100000)/100000;

if ab == 1
    a = 1:M;
    a = a.^(-d);
    a = (a./sum(a)).*R;
end

if ab == 2
    a = 1:M;
    a = exp(-a.*d);
    a = (a./sum(a)).*R;
end

if ab == 3
    a = 1:M;
    a = log10(a+1).^(-d);
    a = (a./sum(a)).*R;
end

if ab == 4
    t = zeros(1,(M+1));
    t(1,1) = -Inf; t(1,(M+1)) = Inf;
    for i = 1:(M-1)
        t(1,(i+1)) = sqrt(2)*erfinv((2/M) + erf(t(1,i)/sqrt(2)));
    end
    
    a = zeros(1,M);
    for i = 1:M
        a(1,i) = (M/sqrt(2*pi)*(exp(-0.5*t(1,i)^2)-exp(-0.5*t(1,(i+1))^2)));
    end
    a = exp(a*d);
    a = (a./sum(a)).*R;
end

H(:,1) = 1 - exp(-a.*((r-o)/L));

for q = 1:length(E)
    E(1,q) = a*(q.*(H.^(q-1)).*((1-H).^2));
    Z(1,q) = a*(q.*(H.^(q-1)).*((1-H).^2).*(1-(q.*(H.^(q-1)).*((1-H).^2))));
end

y = sum(((C-E).^2)./Z);

end

% The Genetic Algorithm iterator defined as below

function [x,fval,ef, runTime] = ga_iteration(i,j,abm,res,LB,UB,n_partitions,csp, cfg)
opt = gaoptimset('TolCon',0,'TolFun',0,'Display','off','PopInitRange',[LB;UB],'UseParallel', true, 'Vectorized', 'off');
IntCon = [1,2,3];

tic;
rng('default');

objective_fn = @(x)cost_fn(x, cfg, csp, abm, res);
[x,fval,ef] = ga(objective_fn,3,[],[],[],[],LB,UB,[],IntCon,opt);
runTime = toc;
x = x.*[res(1), res(2), res(3)];
row_n = (i-1)*n_partitions + j;
fprintf('%i\t%i\t%i\t%i\t%.5f\t%i\t%e\t%i\t%.4f\n',(row_n),j,x(1),x(2),x(3),abm,fval,ef,runTime);
end

function y = calculateEvenness(ab_model, model_parameter,richness)

if (ab_model ==1)
a = 1:richness;
f= a.^(-model_parameter);
f = (f./sum(f));
y = (-1*sum(log(f).*f))/log(richness);
return;
end


if (ab_model==2)
a = 1:richness;
    f=exp(-a.*model_parameter) ;
    f = (f./sum(f));
    y = -1*sum(log(f).*f)/log(richness); 
    return;
end
if (ab_model==3)
a = 1:richness;
    f=log10(a+1).^(-model_parameter);
    f = (f./sum(f));
    y = -1*sum(log(f).*f)/log(richness);
    return;

end
if (ab_model==4)
t = zeros(1,(richness+1));
    t(1,1) = -Inf; t(1,(richness+1)) = Inf;
    for i = 1:(richness-1)
        t(1,(i+1)) = sqrt(2)*erfinv((2/richness) + erf(t(1,i)/sqrt(2)));
    end

    a1 = zeros(1,richness);
    for i = 1:richness
        a1(1,i) = (richness/sqrt(2*pi)*(exp(-0.5*t(1,i)^2)-exp(-0.5*t(1,(i+1))^2)));
    end

    f=exp(a1*model_parameter);
    f = (f./sum(f));
    y = -1*sum(log(f).*f)/log(richness); 
return;
end
end