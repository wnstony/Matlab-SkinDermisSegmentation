% Main file

clear;
close all;
warning off;

addpath(genpath('D:\Matlab\Library\Bouma toolbox'));
addpath(genpath('D:\Matlab\Library\Xinyu toolbox'));
addpath(genpath('D:\Matlab\Library\Colormaps'));
addpath(genpath('D:\Matlab\Library\subaxis'));
addpath(genpath('D:\Matlab\Library\ExportFig'));
addpath(genpath('D:\Matlab\Library\BlandAltman'));
addpath(genpath('d:\Matlab\Library\Symmetrization'));

% Data files
addpath(genpath('..\Data'));

Files{1} = '[PhantomPackage][03-21-2018_15-17-52]';
Files{2} = '[PhantomPackage][03-21-2018_15-18-26]';

fileNum = numel(Files);
caption = cell(fileNum);

dzReal = 10;% in um
fwx = 8;
fwy = 8;
fwz = 8;
NfftBase = 1024;
NfftFold = 1;
Nfft = NfftBase*NfftFold;

% Birefringence calculation related parameters
n = 1.33;
pixDistance = 6.0176/NfftFold/n; % 6.0176um for Nout=1024 (NfftBase). In air.

dz = dzReal/pixDistance;
dz = round((dz-2)/2)*2+2; % round to nearest even integer

dnLimit = 3e-3;
%UADlimit = maxDn*2*pi*dzReal/1.3; % 1.3 um is central wavelength
UADLimit = 2*dnLimit*360/1.3; % Double pass retardation in deg/um. (1.3 um is central wavelength)
%UADlimit = 0.25;  % set a maximum retardation in degree/um
%RetLimit = UADlimit*dzReal*pi/180; % rad/retardation exceeds this value will be set to this value

% result statistical analysis
binsUAD = linspace(0,UADLimit,101);
binsDn = linspace(0,dnLimit,101);
binsOA = linspace(-pi,pi,501);
%% intensity image preview
fIndPreview = 12; % !!!!!!!!!!!!!!!!!!!!!!!!!!!index of file for preview
dispLmt = [14,28]; % gray scale limit for log(int)
% logF = readLogFile(Files{fInd});
% frameNum = logF.numImages;
% stBinned.window = 1;
% stBinned.Nout = Nfft;
% [S1,S2] = recstrTom(Files{fIndPreview},[1,1],stBinned);
% int = tom2Int(S1,S2);
% int = log(int);
% intImg = repmat(mat2gray(int,dispLmt),[1,1,3]);
% 
% figure(100);
% imshow(intImg);

%% Decide the surface and bottom range based on the preview
% =========================================================
% Aline range to cut the return trip of galvo. same for each file.
xRng = [1,1024]; 

% specify z (axial) range to save space and time
zRng = [1;1024]; % 2xN matrix, can be different for each file.

% surface search range, can be different for each file.
sfcRipple = repmat(100,[1,fileNum]);

% a rough thickness of the strip, can be different for each file.
sfcDiff = repmat(270,[1,fileNum]);

% surface control points. 2xnxN matrix, can be different for each file.
sfcCtrlPts = cat(3,repmat([xRng(1),xRng(2); ...
                              1574,705],[1,1,fileNum]) );
%  sfcCtrlPts(:,:,1) = [xRng(1),xRng(2);2038,1374];
%  sfcCtrlPts(:,:,2) = [xRng(1),xRng(2);2038,1278];
%  sfcCtrlPts(:,:,3) = [xRng(1),xRng(2);1886,1150];
%  sfcCtrlPts(:,:,4) = [xRng(1),xRng(2);1670,960];
%  sfcCtrlPts(:,:,5) = [xRng(1),xRng(2);1558,802];
%  
%  sfcCtrlPts(:,:,6) = [xRng(1),xRng(2);1794,877];
%  sfcCtrlPts(:,:,8) = [xRng(1),xRng(2);1690,774];
%  sfcCtrlPts(:,:,9) = [xRng(1),xRng(2);1690,774];
%  sfcCtrlPts(:,:,10) = [xRng(1),xRng(2);1470,518];
%  sfcCtrlPts(:,:,11) = [xRng(1),xRng(2);1690,774];
 
%segMethod = repmat(1,[1,fileNum]); % 1 for peakDetection; 2 for intensHighland
% =========================================================
%% Derived parameters
cutAxialNum = zRng(2,:)-zRng(1,:)+1;
cutAlineNum = xRng(2)-xRng(1)+1;

%% show the image with colored search range
% calculate the rough surface line based on sfcCtrlPts(:,:,fInd);
% sfcCtrlPtsTemp = sfcCtrlPts(:,:,fIndPreview);
% [~,ctrlPtsNum] = size(sfcCtrlPtsTemp);
% sfcLineZ = sfcCtrlPtsTemp(2,1);
% for iPt = 1:ctrlPtsNum-1
%     sfcTemp = linspace(sfcCtrlPtsTemp(2,iPt),sfcCtrlPtsTemp(2,iPt+1),(sfcCtrlPtsTemp(1,iPt+1)-sfcCtrlPtsTemp(1,iPt)+1));
%     sfcLineZ = [sfcLineZ,sfcTemp(2:end)];
% end
% sfcLineZ = round(sfcLineZ);
% sfcLineX = sfcCtrlPtsTemp(1,1):sfcCtrlPtsTemp(1,end);
% % find the 4 boundaries
% colorRange(1,:) = sfcLineZ - sfcRipple(fIndPreview);
% colorRange(2,:) = sfcLineZ + sfcRipple(fIndPreview);
% colorRange(3,:) = sfcLineZ + sfcDiff(fIndPreview) - sfcRipple(fIndPreview);
% colorRange(4,:) = sfcLineZ + sfcDiff(fIndPreview) + sfcRipple(fIndPreview);
% % show upper and lower surface search area in different colors
% for iX = 1:numel(sfcLineX)
%     intImg(colorRange(1,iX):colorRange(2,iX),sfcLineX(iX),1) = 0;
%     intImg(colorRange(3,iX):colorRange(4,iX),sfcLineX(iX),3) = 0;
% end
% % show cut boundary in red line
%     % horizontal line 1
%     intImg(zRng(1,fIndPreview):zRng(1,fIndPreview)+3,xRng(1):xRng(2),1) = 1;
%     intImg(zRng(1,fIndPreview):zRng(1,fIndPreview)+3,xRng(1):xRng(2),2:3) = 0;
%     % horizontal line 2
%     intImg(zRng(2,fIndPreview)-3:zRng(2,fIndPreview),xRng(1):xRng(2),1) = 1;
%     intImg(zRng(2,fIndPreview)-3:zRng(2,fIndPreview),xRng(1):xRng(2),2:3) = 0;
%     % vertical line 1
%     intImg(zRng(1,fIndPreview):zRng(2,fIndPreview),xRng(1):xRng(1)+3,1) = 1;
%     intImg(zRng(1,fIndPreview):zRng(2,fIndPreview),xRng(1):xRng(1)+3,2:3) = 0;
%     % vertical line 2
%     intImg(zRng(1,fIndPreview):zRng(2,fIndPreview),xRng(2)-3:xRng(2),1) = 1;
%     intImg(zRng(1,fIndPreview):zRng(2,fIndPreview),xRng(2)-3:xRng(2),2:3) = 0;
%     
% figure(100);
% imshow(intImg);

%% params
stUnbinned.window = 1; % 'stNoBinning' for segmentation
stUnbinned.Nout = Nfft;
stBinned.window = 5; % 'st' for retardation calculation
stBinned.Nout = Nfft;
stBinned.gpu = false; % max capacity: 4096*1024*9*3; cannot be 8192

% binNum = 2*stBinned.window-1;
% S1w = single(zeros(cutAxialNum,cutAlineNum,binNum,3));
% S2w = S1w;

% midFrameInd = round(frameNum/2);
% frameRange = [midFrameInd-ceil(numel(ny)/2),midFrameInd+ceil(numel(ny)/2)];
% selFrameInd = midFrameInd - frameRange(1) - 1;
frameInd = 1;
%% main loop: File & Frame
recDn = zeros(cutAxialNum,cutAlineNum,fileNum);
recRet = zeros(cutAxialNum,cutAlineNum,fileNum);
recOA = zeros(cutAxialNum,cutAlineNum,fileNum);
recDOP = zeros(cutAxialNum,cutAlineNum,fileNum);

fIndOffset = 0;
for fInd = 2 % file loop
    fInd
    % read file notes
    logF = readLogFile(Files{fInd});
    caption(fInd) = cellstr(logF.msmtNotes);
    
    %% birefringence reconstruction    
    % reconstruct binned Stokes and cut ROI
    [S1_Orig,S2_Orig] = recstrTom(Files{fInd+fIndOffset},frameInd,stBinned); % process one frame for each loop
    S1w = S1_Orig(zRng(1):zRng(2),xRng(1):xRng(2),:,:);
    S2w = S2_Orig(zRng(1):zRng(2),xRng(1):xRng(2),:,:);
    clear S1_Orig S2_Orig;
    
    pstruct.fwx = fwx; % x average happens here
    out = ProcessForSymmetry(S1w,S2w,pstruct);
    
    %         StokesIn.Sto1 = S1w; %1024*1024*9*3
    %         StokesIn.Sto2 = S2w;
    
    StokesIn.Sto1 = permute(out.MM(1:3,:,:,:),[2,3,4,1]); %1024*1024*9*3
    StokesIn.Sto2 = permute(out.MM(4:6,:,:,:),[2,3,4,1]);
    params.dz = dz;
    %params.fwz = 0; % Depth average happens here
    frameOut = SpectralBinning_StoIn(StokesIn,out.dop,params);
    
    %% find the surface and bottom line
    [S1,S2] = recstrTom(Files{fInd+fIndOffset},frameInd,stUnbinned);
    int = log(tom2Int(S1,S2));  % following segmentation all based on 'int'
    int = int(zRng(1):zRng(2),xRng(1):xRng(2));
    intImg = repmat(mat2gray(int,dispLmt),[1,1,3]);
    
%     estSfcLine = f_findSurfaceUsingDOP(out.dop,0.7);
%     estSfcLine = reshape(estSfcLine,[1,numel(estSfcLine)]); % make sure it's 1xN matrix
    
    % -- cut sfcCtrlPts, Row1-xRng(1)+1; Row2-zRng(1)+1
    
%     sfcCtrlPtsCut = [1,size(out.dop,2);1,1];
%     if segMethod(fInd) == 1
%         [segmentMask,unifiedThickness] = f_seg_peakDetection(int,dispLmt,estSfcLine,sfcDiff(fInd),sfcRipple(fInd));
%     else
%         [segmentMask,unifiedThickness] = f_seg_intensityHighland(int,dispLmt,estSfcLine,sfcDiff(fInd),sfcRipple(fInd));
%     end
%     
%     maskAlineNum = size(out.dop,2);
    %maskAlineNum = sfcCtrlPtsCut(1,end)-sfcCtrlPtsCut(1,1)+1;
    %% Ret/DOP/OA array in the segmented area
    angleArray = frameOut.ret*dz;
    retArray = angleArray*180/pi/(dz*pixDistance); % double pass ret deg/um
    dnArray = angleArray*1.3/(4*pi)/(dz*pixDistance);
    oaArray = atan(frameOut.PA(:,:,2)./frameOut.PA(:,:,1));
    dopArray = out.dop;
    
%     angleArray = zeros(unifiedThickness,maskAlineNum);
%     retArray = zeros(unifiedThickness,maskAlineNum);
%     dnArray = zeros(unifiedThickness,maskAlineNum);
%     dopArray = zeros(unifiedThickness,maskAlineNum);
%     oaArray = zeros(unifiedThickness,maskAlineNum);
%     for iX = 1:maskAlineNum
%         iXInMask = iX+sfcCtrlPtsCut(1,1)-1;
%         angleArray(:,iX) = rotAngle(segmentMask(:,iXInMask)==1,iXInMask);
%         retArray(:,iX) = ret(segmentMask(:,iXInMask)==1,iXInMask);
%         dnArray(:,iX) = dn(segmentMask(:,iXInMask)==1,iXInMask);
%         oaArray(:,iX) = oa(segmentMask(:,iXInMask)==1,iXInMask);
%         dopArray(:,iX) = out.dop(segmentMask(:,iXInMask)==1,iXInMask);
%     end
%     moreShrink1 = 70;
%     moreShrink2 = 70;
%     moreShrink = moreShrink1 + moreShrink2;
%     
%     angleArray = angleArray(moreShrink1+1:end-moreShrink2,:);
%     retArray = retArray(moreShrink1+1:end-moreShrink2,:);
%     dnArray = dnArray(moreShrink1+1:end-moreShrink2,:);
%     oaArray = oaArray(moreShrink1+1:end-moreShrink2,:);
%     dopArray = dopArray(moreShrink1+1:end-moreShrink2,:);

    %dnArray = retArray*1.3/360/2;    
    retArray(retArray>UADLimit) = UADLimit;
    dnArray(dnArray>dnLimit) = dnLimit;
    oaArray(isnan(oaArray)) = 0;
    oaArray = real(oaArray);
    
    %% recording for each file
%     recAlineNum(fInd) = maskAlineNum;
%     recThickness(fInd) = unifiedThickness;
    recDn(:,:,fInd) = dnArray;
    recRet(:,:,fInd) = retArray;
    recOA(:,:,fInd) = oaArray;
    recDOP(:,:,fInd) = dopArray;
    %% Statistic
%     histDn(:,fInd) = histc(reshape(dnArray,[numel(dnArray),1]),binsDn);
%     histUAD(:,fInd) = histc(reshape(retArray,[numel(retArray),1]),binsUAD);
%     histOA(:,fInd) = histc(reshape(oaArray,[numel(oaArray),1]),binsOA);
    
    %% Draw results
    spacing = 0.04;
    padding = 0.01;
    mgSide = 0.03;
    mgTop = 0.08;
    mgBtm = 0.03;
    yNum = 2;
    xNum = 2;
    f = figure(1);
    set(f,'Position', [100 50 1900 1250])
    
%     subaxis(yNum,xNum,1,'Spacing',spacing,'Padding',padding,'Margin',margin);
%     imagesc(ret,[0,UADlimit])
%     colormap(flipud(ametrine))
%     colorbar;
    
    ax1 = subaxis(yNum,xNum,1,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
    imshow(intImg);
    title('Intensity');
    set(gca,'fontsize',15);
    
    ax2 = subaxis(yNum,xNum,2,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
    imagesc(dnArray,[0,dnLimit])
    colormap(flipud(ametrine))
    colorbar;
    title('retardation (delta n)');
    set(gca,'fontsize',15);
    
    ax3 = subaxis(yNum,xNum,3,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
    imagesc(dopArray,[0,1])
    colormap(flipud(ametrine))
    colorbar;
    title('DOP');
    set(gca,'fontsize',15);
    
    ax4 = subaxis(yNum,xNum,4,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
    imagesc(oaArray,[-pi,pi])
    colormap(flipud(ametrine))
    colorbar;
    title('OA');
    set(gca,'fontsize',15);
   
    
    suptitle(caption(fInd));
    set(gcf,'color','w');
    
    imgName = sprintf(char(strcat('%d_',caption(fInd),'.fig')),fInd);
    imgName = convertCharsToStrings(imgName);
%     export_fig imgName -q101
    savefig(imgName)
    
    %     end % END of frame loop
end % END of file loop

% f = figure(2);
% set(f,'Position', [100 50 1900 1250])
% plot(binsDn,histDn(:,1:7))
% xlim(ax2,[0 dnLimit]);
% lgd = legend(caption(1:7));
% lgd.FontSize = 26;
% 
% savefig('unannealed stretched compare.fig')
% 
% f = figure(3); % unstretched compare
% set(f,'Position', [100 50 1900 1250])
% ax1 = subaxis(yNum,xNum,1,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
% imagesc(recDn(1:recThickness(8),:,8),[0,dnLimit])
% colormap(flipud(ametrine))
% colorbar;
% title('unannealed retardation (deg/um)');
% set(gca,'fontsize',15);
% ax2 = subaxis(yNum,xNum,2,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
% imagesc(recDn(1:recThickness(9),:,9),[0,dnLimit])
% colormap(flipud(ametrine))
% colorbar;
% title('annealed retardation (deg/um)');
% set(gca,'fontsize',15);
% ax3 = subaxis(yNum,xNum,3,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',mgTop);
% plot(binsDn,histDn(:,8))
% hold on;
% plot(binsDn,histDn(:,9))
% xlim(ax3,[0 dnLimit]);
% title('retardation distribution');
% set(gca,'fontsize',15);
% lgd = legend({'Unannealed','annealed'});
% lgd.FontSize = 26;
% ax4 = subaxis(yNum,xNum,4,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',mgTop);
% plot(squeeze(mean(recDn(1:recThickness(8),:,8:9),2)))
% ylim(ax4,[0 dnLimit]);
% title('retardation mean across depth');
% set(gca,'fontsize',15);
% lgd = legend({'Unannealed','annealed'});
% lgd.FontSize = 26;
% 
% suptitle('unstretched compare');
% set(gcf,'color','w');
% 
% savefig('unstretched compare.fig')
% 
% yNum = 1;
% f = figure(4); % heating compare
% set(f,'Position', [50 50 2050 1050])
%     ax1 = subaxis(yNum,xNum,1,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
%         plot(binsDn,histDn(:,1),'LineWidth',1.6)
%         hold on;
%         plot(binsDn,histDn(:,2),'LineWidth',1.6)
%         plot(binsDn,histDn(:,3),'LineWidth',1.6)
%         plot(binsDn,histDn(:,7),'LineWidth',1.6)
%         xlim(ax1,[0 dnLimit]);
%         title('retardation distribution');
%         set(gca,'fontsize',15,'linewidth',2);
%     ax2 = subaxis(yNum,xNum,2,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
%         plot(squeeze(mean(recDn(1:recThickness(1),:,1),2)),'LineWidth',1.6)
%         hold on;
%         plot(squeeze(mean(recDn(2:recThickness(2),:,2),2)),'LineWidth',1.6)
%         plot(squeeze(mean(recDn(3:recThickness(3),:,3),2)),'LineWidth',1.6)
%         x = 11:194;
%         plot(x,squeeze(mean(recDn(7:recThickness(7),:,7),2)),'LineWidth',1.6)
%         ylim(ax2,[0 dnLimit]);
%         title('retardation mean across depth');
%         set(gca,'fontsize',15,'linewidth',2);
%         lgd = legend([caption(1:3),caption(7)]);
%         lgd.FontSize = 26;
% suptitle('heating compare','fontsize',25);
% set(gcf,'color','w');
% savefig('heating compare.fig')
% 
% yNum = 1;
% f = figure(5); % cooling compare
% set(f,'Position', [50 50 2050 1050])
%     ax1 = subaxis(yNum,xNum,1,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
%         plot(binsDn,histDn(:,1),'LineWidth',1.6)
%         hold on;
%         plot(binsDn,histDn(:,4),'LineWidth',1.6)
%         plot(binsDn,histDn(:,5),'LineWidth',1.6)
%         plot(binsDn,histDn(:,6),'LineWidth',1.6)
%         xlim(ax1,[0 dnLimit]);
%         title('retardation distribution');
%         set(gca,'fontsize',15,'linewidth',2);
%     ax2 = subaxis(yNum,xNum,2,'MarginTop',mgTop,'MarginBottom',mgBtm,'MarginLeft',mgSide,'MarginRight',mgSide,'Spacing',spacing,'PaddingTop',padding);
%         plot(squeeze(mean(recDn(1:recThickness(1),:,1),2)),'LineWidth',1.6)
%         hold on;
%         plot(squeeze(mean(recDn(4:recThickness(4),:,4),2)),'LineWidth',1.6)
%         plot(squeeze(mean(recDn(5:recThickness(5),:,5),2)),'LineWidth',1.6)
%         x = 11:194;
%         plot(squeeze(mean(recDn(6:recThickness(6),:,6),2)),'LineWidth',1.6)
%         ylim(ax2,[0 dnLimit]);
%         title('retardation mean across depth');
%         set(gca,'fontsize',15,'linewidth',2);
%         lgd = legend([caption(1),caption(4:6)]);
%         lgd.FontSize = 26;
% suptitle('Cooling compare','fontsize',25);
% set(gcf,'color','w');
% savefig('Cooling compare.fig')
% 
% %% statistic params
% for i = 1:9
%     temp = squeeze(recDn(1:recThickness(i),:,i));
%     temp = reshape(temp,[numel(temp),1]);
%     stt(i,1) = mean(temp(:));
%     stt(i,2) = var(temp(:));
% end




