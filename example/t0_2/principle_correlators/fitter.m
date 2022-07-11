%%%%%% read in principle correlator %%%%%%%%%%%%%%%
prin_corr_filename = 'principle_correlator_6.jack';
svd_cut = 1e-12;
t0 =2;
tmin =3;
tmax = 20;
formatSpec = '%f';
addpath("/home/chris/projects/GEVP/fitters/");
file_ = fopen(prin_corr_filename,'r');
format long;
prin_corr = fscanf(file_,formatSpec);

%%%%%%%%%%% Need to reformat principle_correlator into matrix %%%%%

Ncfgs = prin_corr(1); %%% Number of configurations
Nt = prin_corr(2);    %%% Number of timeslices

ts = prin_corr(6:2:(2*Nt+4));
vals = prin_corr(7:2:end);
vals = reshape(vals, Nt, Ncfgs);

%%%%%%%%%%% Build covariance matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cov = zeros(Nt,Nt);
average = mean(vals,2);
for t1 = 1:Nt

    for t2 = 1:Nt
        cov_little =0;
        for cfg = 1:Ncfgs
            cov_little = cov_little + (((vals(t1,cfg)-average(t1))*(vals(t2,cfg)-average(t2))));
        end
        cov_little = cov_little/(Ncfgs);
        cov_little = cov_little/(Ncfgs-1);
        Cov(t1,t2) = cov_little;
        
    end
end
InvCov = pinv(Cov,svd_cut);

%%%%%%%%%%% Create jacknife average ensemble %%%%%%%%%%%%%%%%%%%%
average = average.*ones(Nt,Ncfgs);
vals_jack = (Ncfgs/(Ncfgs-1))*average - (1/(Ncfgs-1))*vals;


%%%%%%%%%%% Minimize Chisquare for each average jackknife configuration

options = optimset('MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-4);
a = zeros(3,Ncfgs);
for cfg = 1:Ncfgs
    a(:,cfg)=fminsearch(@(x) chisq_function(tmin,tmax,t0,InvCov,vals_jack(:,cfg),x), [[0.2,0.3,1]],options);
end

fit_params_average = mean(a,2);

%%%%%%%%%%% undo jacknife bias %%%%%%%%%%%%%%%%%%%%%%%%%


fit_params = Ncfgs*fit_params_average.*ones(3,Ncfgs)-(Ncfgs-1)*a;

errors = sqrt(var(fit_params,0,2)/(Ncfgs-1));
fileID = fopen('mass.jack','w');
formatSpec= ' 0 %6.6f\n';
fprintf(fileID,'%d 1 0 0 1\n',Ncfgs);
fprintf(fileID,formatSpec,fit_params(2,:));
fclose(fileID);

chisq_per_dof = chisq_function(tmin,tmax,t0,InvCov,mean(vals,2), fit_params_average);

times = [0:Nt-1];

prin_corr_errs = sqrt(var(vals,0,2));
prin_corr_central = mean(vals,2);
prin_corr_central(t0+1) = 1;
prin_corr_errs(t0+1)=0;

prin_corr_weight = zeros(Nt,1); prin_corr_errs_weight = zeros(Nt,1);
for t = 0:Nt-1
    prin_corr_weight(t+1) = prin_corr_central(t+1)*exp(fit_params_average(2)*(t-t0));
    prin_corr_errs_weight(t+1) = prin_corr_errs(t+1)*exp(fit_params_average(2)*(t-t0));
end

%%%%%%%% plot data %%%%%%%%%%%%%%
e = errorbar(times,prin_corr_weight,prin_corr_errs_weight);
e.LineStyle = 'none';
e.Marker = 'o';
e.MarkerSize = 3;
e.MarkerFaceColor = 'w';
axis([0 tmax+3 0 4]);
set(gcf, 'renderer', 'painters')

hold on
%%%%%%% plot fit before fit range %%%%%%%%%%%
tfirst = linspace(0,tmin);
central = exp(fit_params_average(2).*(tfirst-t0)).*((1-fit_params_average(1)).*exp(-fit_params_average(2).*(tfirst-t0))+fit_params_average(1).*exp(-fit_params_average(3).*(tfirst-t0)));
fc = plot(tfirst,central,'Color','b');
set(gcf, 'renderer', 'painters')

hold on

plus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(errors(1)+fit_params_average(1))).*exp(-(errors(2)+fit_params_average(2)).*(tfirst-t0))+(errors(1)+fit_params_average(1)).*exp(-(errors(3)+fit_params_average(3)).*(tfirst-t0)));
fp = plot(tfirst,plus,'Color','b');
set(gcf, 'renderer', 'painters')

hold on

minus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(-errors(1)+fit_params_average(1))).*exp(-(-errors(2)+fit_params_average(2)).*(tfirst-t0))+(-errors(1)+fit_params_average(1)).*exp(-(-errors(3)+fit_params_average(3)).*(tfirst-t0)));
fm = plot(tfirst,minus,'Color','b');
set(gcf, 'renderer', 'painters')

hold on

%%%%%%%%% plot fit %%%%%%%%%%%%%%%%%

tfirst = linspace(tmin,tmax);
central = exp(fit_params_average(2).*(tfirst-t0)).*((1-fit_params_average(1)).*exp(-fit_params_average(2).*(tfirst-t0))+fit_params_average(1).*exp(-fit_params_average(3).*(tfirst-t0)));
fc = plot(tfirst,central,'Color','r');
set(gcf, 'renderer', 'painters')

hold on

plus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(errors(1)+fit_params_average(1))).*exp(-(errors(2)+fit_params_average(2)).*(tfirst-t0))+(errors(1)+fit_params_average(1)).*exp(-(errors(3)+fit_params_average(3)).*(tfirst-t0)));
fp = plot(tfirst,plus,'Color','r');
set(gcf, 'renderer', 'painters')

hold on

minus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(-errors(1)+fit_params_average(1))).*exp(-(-errors(2)+fit_params_average(2)).*(tfirst-t0))+(-errors(1)+fit_params_average(1)).*exp(-(-errors(3)+fit_params_average(3)).*(tfirst-t0)));
fm = plot(tfirst,minus,'Color','r');
set(gcf, 'renderer', 'painters')

hold on

%%%%%%%%% plot fit after fit range %%%%%%%%%%%%%%%%%

tfirst = linspace(tmax,Nt);
central = exp(fit_params_average(2).*(tfirst-t0)).*((1-fit_params_average(1)).*exp(-fit_params_average(2).*(tfirst-t0))+fit_params_average(1).*exp(-fit_params_average(3).*(tfirst-t0)));
fc = plot(tfirst,central,'Color','b');
set(gcf, 'renderer', 'painters')

hold on

plus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(errors(1)+fit_params_average(1))).*exp(-(errors(2)+fit_params_average(2)).*(tfirst-t0))+(errors(1)+fit_params_average(1)).*exp(-(errors(3)+fit_params_average(3)).*(tfirst-t0)));
fp = plot(tfirst,plus,'Color','b');
set(gcf, 'renderer', 'painters')

hold on

minus = exp(fit_params_average(2).*(tfirst-t0)).*((1-(-errors(1)+fit_params_average(1))).*exp(-(-errors(2)+fit_params_average(2)).*(tfirst-t0))+(-errors(1)+fit_params_average(1)).*exp(-(-errors(3)+fit_params_average(3)).*(tfirst-t0)));
fm = plot(tfirst,minus,'Color','b');
set(gcf, 'renderer', 'painters')

hold on
txt1 = ['chisq per dof = ' num2str(chisq_per_dof)];
txt2 = ['Mass = ' num2str(fit_params_average(2)) '\pm' num2str(errors(2))];
text(tmin,3,txt1)
text(tmin,2.5,txt2)

%%%%%%%%%% Function declerations %%%%%%%%%%%%%%%%%%%

function correlator = lambda_corr(tmin,tmax,t0,A,m0,m1)
  
   correlator = zeros(tmax-tmin,1);
   for t = 0:(tmax-tmin)
       correlator(t+1) = (1-A)*exp(-m0*(t+tmin-t0))+A*exp(-m1*(t+tmin-t0));
   end
end

function chisq = chisq_function(tmin,tmax,t0,InvCov,jack,y)
%function chisq = chisq_function(x)
%%%%%%%% x[0]=A, x[1]=m0, x[2] = m1
   
   chisq = 0;
   lambda = lambda_corr(tmin,tmax,t0,y(1),y(2),y(3));
   for t1=0:tmax-tmin

       for t2 = 0:tmax-tmin
          
           chisq = chisq + (lambda(t1+1)-jack(tmin+t1+1))*InvCov(t1+tmin+1,t2+tmin+1)*(lambda(t2+1)-jack(tmin+t2+1));
       end
   end
   chisq = chisq /(tmax-tmin-3);

end