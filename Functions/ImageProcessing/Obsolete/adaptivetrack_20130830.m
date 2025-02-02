function [tracked,newdaughters,nuc_label]=adaptivetrack(prevdata,curdata,nuc_raw,nuc_label,nucr,extractmask,jitx,jity)
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
bordermask=zeros(size(nuc_mask));
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=6*nucr;         %farthest jump (pixels)
mitosisrad=3*nucr; %farthest jump to daughter (pixels)
xprev=prevdata(:,1); yprev=prevdata(:,2); areaprev=prevdata(:,3); massprev=prevdata(:,4);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcur=numel(xcur);
mergeid=zeros(numcur,1);
borderflag=zeros(numcur,1);
%%% detect merges in current frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numcur
	neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    %%% find closest neighbor and check mass difference %%%%%%%%%%%%%%%%%%
    [~,cidx]=min(sqrt((xprev(neighbors)-xcur(i)).^2+(yprev(neighbors)-ycur(i)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(i)-massprev(match))/massprev(match);
    if massdiff>0.2 %possible merge detected
        mergeid(i)=1;
        %%% attempt to segment deflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [r,c]=find(nuc_label_mask==i);
        coorset=[c,r];  %adjust to x-y convention
        [order,status]=orderperimeter([c,r]);
        if status==0    %unable to order perimeter
            fprintf('unable to order perimeter\n');
            continue;
        end
        orderedset=coorset(order,:);
        [bordermask,borderflag(i)]=splitdeflections(orderedset,bordermask,nucr);
    end
end
%%% assign un-merged cells new IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitid=find(mergeid & borderflag);
mergemaskorg=ismember(nuc_label,splitid);
mergemask=mergemaskorg & ~bordermask;
mergemask=~bwmorph(~mergemask,'diag');
newlabels=bwlabel(mergemask);
%%% extract features of new cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_info=struct2cell(regionprops(newlabels,nuc_raw,'Area','Centroid','MeanIntensity')');
new_area=squeeze(cell2mat(new_info(1,1,:)));
new_center=squeeze(cell2mat(new_info(2,1,:)))';
new_density=squeeze(cell2mat(new_info(3,1,:)));
new_bg=getbackground_H2B(nuc_label,new_center,nuc_raw,nucr,0.25);
new_density=new_density-new_bg;
new_mass=new_density.*new_area;
%%% update all data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xcur=[xcur;new_center(:,1)]; xcur(splitid)=NaN;
ycur=[ycur;new_center(:,2)]; ycur(splitid)=NaN;
areacur=[areacur;new_area];  areacur(splitid)=NaN;
masscur=[masscur;new_mass];  masscur(splitid)=NaN;
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%
numprevorg=numel(xprev);
previd=find(~isnan(xprev));
numprev=numel(previd);
newdaughters=ones(numcur,1)*NaN;
xprev=xprev(previd); yprev=yprev(previd); massprev=massprev(previd);
prevmatch=ones(numprev,1)*NaN;
%distances=zeros(numprev,1);
%masschanges=zeros(numprev,1);
nextnearestflag=zeros(numprev,1);
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    %%% in case of zero neighbors: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    %%% in case of only one neighbor: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(neighbors)==1
        if abs(massdiff)<0.2
            prevmatch(i)=neighbors;
            %distances(i)=dist;
            %masschanges(i)=massdiff;
        end
        continue;
    end
    %%% for multiple neighors, find closest 2 neighbors %%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    bestidx=didx(1);
    nextidx=didx(2);
    candidateidx=[bestidx;nextidx];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx);
    %%%%%% check for daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    daughtercheck=massdiff>-0.55 & massdiff<-0.40 & dist<mitosisrad;
    if sum(daughtercheck)==2
        newdaughters(candidates)=previd(i);
        prevmatch=[prevmatch;candidates]; %only used to resolve conflicts
        %distances=[distances;0;0];        %daughter call always wins
        %masschanges=[masschanges;0;0];
        continue;
    end
    %%%%%% disqualify inconsistent candidates %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    consistencycheck=abs(massdiff)<0.2;
    if consistencycheck(1)==1
        prevmatch(i)=candidates(1);
        %distances(i)=dist(1);
        %masschanges(i)=massdiff(1);
    elseif consistencycheck(2)==1
        nextnearestflag(i)=1;
        prevmatch(i)=candidates(2);
        %distances(i)=dist(2);
        %masschanges(i)=massdiff(2);
    end
end
%%% resolve conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempsort=sort(prevmatch);
curconflicts=unique(tempsort(logical([diff(tempsort)==0;0])));
for i=curconflicts'
    prevcells=find(prevmatch==i);
    flagcheck=nextnearestflag(prevcells);
    goodprevcells=prevcells(flagcheck==0);
    if numel(goodprevcells)==1
        badprevcells=prevcells(flagcheck==1);
        prevmatch(badprevcells)=NaN;
    else
        conflictdaughter=prevcells(prevcells>numprev);
        if numel(conflictdaughter)==1 %single daughter should win
            prevcells(prevcells~=conflictdaughter)=NaN;
        elseif numel(conflictdaughter)>1 %must remove the newdaughters assignment
            newdaughters(conflictdaughter)=NaN;
            prevmatch(prevcells)=NaN;
        else
            prevmatch(prevcells)=NaN;
        end
    end
end
%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prevmatch=prevmatch(1:numprev); %remove daughters
tracked=ones(numprevorg,4)*NaN;
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempidx=find(~isnan(prevmatch));
matchidxprev=previd(tempidx);
matchidxcur=prevmatch(tempidx);
tracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nontrackedidx=find(~ismember((1:numcur)',prevmatch));
nontracked=[xcur(nontrackedidx) ycur(nontrackedidx) areacur(nontrackedidx) masscur(nontrackedidx)];
mothers=newdaughters(~isnan(newdaughters));
goodnewdaughters=find(~isnan(newdaughters));
newdaughteridx=ismember(nontrackedidx,goodnewdaughters);
tracked=[tracked;nontracked];
newdaughters=ones(numel(nontrackedidx),1)*NaN;
newdaughters(newdaughteridx)=mothers;
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relabelidx=ones(numcur,1)*NaN;
relabelidx(matchidxcur)=matchidxprev;
relabelidx(nontrackedidx)=numprevorg+1:numprevorg+numel(nontrackedidx);
nuc_label(mergemaskorg>0)=0;
newlabels=newlabels+numcur;       %offset so that all label IDs are new
newlabels(newlabels==numcur)=0;
nuc_label=nuc_label+newlabels;
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
[height,width]=size(nuc_mask);
dxminprev=max([round(xprev(i)-winrad) 1]); dxmaxprev=min([round(xprev(i)+winrad) width]);
dyminprev=max([round(yprev(i)-winrad) 1]); dymaxprev=min([round(yprev(i)+winrad) height]);
dxmincur=round(dxminprev-jitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-jitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-jity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-jity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=mat2gray(dbmaskprev);
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(dbimage);
%}