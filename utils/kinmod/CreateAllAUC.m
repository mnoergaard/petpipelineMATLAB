%
files=dir('*.tac');
for i=1:length(files)
	fprintf('Processing (%i): %s\n',i,files(i).name);
	CreateAUC(files(i).name);
end
