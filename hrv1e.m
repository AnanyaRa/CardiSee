clear;
close all;
indata=readtable('data.dat');
X=[indata.cm20 indata.z3e20 indata.SDNN indata.PIP indata.LFHFr];
Y=ordinal(indata.Class);

%-Vizualization
classmat=table2array(indata(1:end,2));
cm20pmat=table2array(indata(1:end,3));
z3e20mat=table2array(indata(1:end,4));
sdnnmat=table2array(indata(1:end,5));
pipmat=table2array(indata(1:end,6));
lfhfrmat=table2array(indata(1:end,7));

figure
gscatter(pipmat,lfhfrmat,classmat,'rgb','osd')
grid on;
%-
figure
gscatter(sdnnmat,cm20pmat,classmat,'rgb','osd');
set(gca, 'YScale', 'log');
grid on;

%- gplotmatrix
figure
xvar=[cm20pmat z3e20mat sdnnmat];
yvar=[pipmat lfhfrmat];
gplotmatrix(xvar,yvar,classmat)
[~,ax]=gplotmatrix(xvar,yvar,classmat);
xlabel(ax(1,1),'CM20'); xlabel(ax(1,2),'Z3e20'); xlabel(ax(1,3),'SDNN');
xlabel(ax(2,1),'CM20'); xlabel(ax(2,2),'Z3e20'); xlabel(ax(2,3),'SDNN');
ylabel(ax(1,1),'PIP'); ylabel(ax(2,1),'LF/HF ratio');

for i=1:2
    for j=1:2
    set(ax(1,j),'xscale','log');
    end
end

%- 3D scatter plot
[~,~,id] = unique(classmat);
colors = 'rgb';
markers = 'osd';
figure
for idx = 1 : 3
    data = X(id == idx,:);
    plot3(data(:,1), data(:,2), data(:,3), [colors(idx) markers(idx)]);
    hold on;
end
set(gca, 'XScale', 'log');
title('Scatter Plot: cm20p vs Z3e vs SDNN')
xlabel('cm20p'); ylabel('Z3e'); zlabel('SDNN');
legend({'x=CM','y=Z3e','z=SDNN'});
grid;
%--
figure

for idx = 1 : 3
    data = X(id == idx,:);
    plot3(data(:,1), data(:,3), data(:,5), [colors(idx) markers(idx)]);
    hold on;
end

set(gca, 'XScale', 'log');
title('Scatter Plot: cm20p vs SDNN vs. LFHF ratio');
xlabel('cm20p'); ylabel('SDNN'); zlabel('LFHF ratio');
legend({'x=CM','y=SDNN','z=LFHF ratio'});
grid;
%--
figure

for idx = 1 : 3
    data = X(id == idx,:);
    plot3(data(:,3), data(:,4), data(:,5), [colors(idx) markers(idx)]);
    hold on;
end

title('Scatter Plot: SDNN vs PIP vs LFHF rato');
xlabel('SDNN'); ylabel('PIP'); zlabel('LFHF ratio');
legend({'x=SDNN','y=PIP','z=LFHF ratio'});
grid;

%- Extract data by one class
indata.Class=categorical(indata.Class);
grpstats(indata(:,2:end),'Class',@median)
%--

arrdata=indata(indata.Class=='ARR',:);
% boxplot(table2array(arrdata(:,3:7)))
% vs= HRV.violinplot(table2array(arrdata(:,3:7)))

figure
subplot(3,1,1)
heatmap(corr(table2array(arrdata(:,3:end))),'XData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"],'YData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"]);
axp=struct(gca);
axp.Axes.XAxisLocation = 'top';
caddata=indata(indata.Class=='CAD',:);

subplot(3,1,2)
heatmap(corr(table2array(caddata(:,3:end))),'XData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"],'YData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"]);
axp=struct(gca);
axp.Axes.XAxisLocation = 'top';
nordata=indata(indata.Class=='NOR',:);

subplot(3,1,3)
heatmap(corr(table2array(nordata(:,3:end))),'XData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"],'YData',["cm20p" "z3e20" "SDNN" "PIP" "LFHFr"]);
axp=struct(gca);
axp.Axes.XAxisLocation = 'top';

subplot(3,1,1)
title('Class ARR correlations')
subplot(3,1,2)
title('Class CAD correlations')
subplot(3,1,3)
title('Class NOR correlations')

%-
%- Correlation Coefficient Analysis of raw data
indata2=table2array(indata(:,3:end));
%z1=zscore(indata2,[ ],1);
ccoeff=corrcoef(indata2);
figure
%heatmap(ccoeff)
heatmap(corr(indata2))

%-
% https://www.mathworks.com/help/stats/examples/credit-rating-by-bagging-decision-trees.html
%--
leaf = [1 5 10];
nTrees = 50;
rng(9876,'twister');
savedRng = rng; % save the current RNG settings
figure
color = 'bgry';

for ii = 1:length(leaf)
   % Reinitialize the random number generator, so that the
   % random samples are the same for each leaf size
   rng(savedRng);
   % Create a bagged decision tree for each leaf size and plot out-of-bag
   % error 'oobError'
   %b = TreeBagger(nTrees,X,Y,'OOBPred','on',...
   %                          'CategoricalPredictors',5,...
   %                          'MinLeaf',leaf(ii));
   b = TreeBagger(nTrees,X,Y,'OOBPred','on',...
                             'MinLeaf',leaf(ii));
   plot(b.oobError,color(ii));
   hold on;
end

xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');
legend({'1', '5','10'},'Location','NorthEast');
title('Classification Error for Different Leaf Sizes');
grid on;
hold off;

%nTrees = 30;
nTrees = 30;
leaf = 1;
rng(savedRng);
%b = TreeBagger(nTrees,X,Y,'OOBVarImp','on','CategoricalPredictors',5,'MinLeaf',leaf);
 b = TreeBagger(nTrees,X,Y,'OOBVarImp','on','MinLeaf',leaf);                     
view(b.Trees{1},'Mode','graph')

figure
b_err=b.OOBPermutedVarDeltaError;
%bar(b.OOBPermutedVarDeltaError);
bar(b_err);
b_err_names={'cm20p' 'z3e20' 'SDNN' 'PIP' 'LFHFr'};
set(gca,'xticklabel',b_err_names);
xlabel('Feature number');
ylabel('Out-of-bag feature importance');
title('Feature importance results');
grid on;
hold off;

oobErrorFullX = b.oobError;
b = b.compact;

newDS=readtable('data_1.dat');
%newDS=readtable('data.dat');
[predClass,classifScore] = b.predict([newDS.cm20 newDS.z3e20 newDS.SDNN newDS.PIP newDS.LFHFr]);

for i = 1:size(newDS,1)
   fprintf('Patient %s:\n',char(newDS.Dataset(i)));
   fprintf('   cm20    = %5.2f\n',newDS.cm20(i));
   fprintf('   z3e    = %5.2f\n',newDS.z3e20(i));
   fprintf('   PIP = %5.2f\n',newDS.PIP(i));
   fprintf('   LFHFr = %5.2f\n',newDS.LFHFr(i));
   fprintf('   SDNN = %5.2f\n',newDS.SDNN(i));
   fprintf('   Actual Cat : %s\n',char(newDS.Class(i)));
   fprintf('   Predicted Cat : %s\n',predClass{i});
   %fprintf('   Classification score : \n');
   %for j = 1:length(b.ClassNames)
   %   if (classifScore(i,j)>0)
   %      fprintf('      %s : %5.4f \n',b.ClassNames{j},classifScore(i,j));
   %   end
   %end
end

%https://www.mathworks.com/help/stats/classification-treeBagger-examples.html#br0g6t1-1
%b5v = TreeBagger(nTrees,X,Y,'OOBVarImp','on','OOBPred','on');
%gPosition = find(strcmp('NOR',b5v.ClassNames));
%[Yfit,Sfit] = oobPredict(b5v);
%[fpr,tpr,~,auc] = perfcurve(b5v.Y,Sfit(:,gPosition),'NOR');
%figure
%plot(fpr,tpr)
%grid on;
%xlabel('False Positive Rate')
%ylabel('True Positive Rate')
%text(0.5,0.25,strcat('AUC=',num2str(auc)),'EdgeColor','k');
%title('ROC curve NOR, predicted vs. actual rating');

%[predClass,classifScore] = b.predict([indata.cm20 indata.z3e20 indata.SDNN indata.PIP indata.LFHFr])

classnames = b.ClassNames;
predDS = [table(newDS.Dataset,predClass),array2table(classifScore)];

predDS.Properties.VariableNames = {'Dataset','PredRating',classnames{:}};
writetable(predDS,'PredictedRatings.dat');
predDS = readtable('PredictedRatings.dat');

C = confusionmat(predDS.PredRating,newDS.Class,'order',{'ARR' 'CAD' 'NOR'})
Cperc=diag(sum(C,2))\C

% Sensitivity and Specifivity
cmat=[C(1,1)+C(2,2) C(2,3);C(3,2) C(3,3)];
selec=cmat(1,1)/sum(cmat(1,:));
spec=cmat(2,2)/sum(cmat(2,:));

% Alternate method for AUC
[xVal,yVal,~,auc] = perfcurve(newDS.Class,predDS.NOR,'NOR');
figure
plot(xVal,yVal);
grid on;
xlabel('False positive rate');
ylabel('True positive rate');
text(0.5,0.25,strcat('AUC=',num2str(auc)),'EdgeColor','k');
title('ROC curve NOR, predicted vs. actual rating');

