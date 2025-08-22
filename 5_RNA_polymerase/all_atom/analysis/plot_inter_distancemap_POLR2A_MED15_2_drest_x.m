clear all;

constructs={'POLR2A_IDR1_15_repeats_w_MED15_IDR2_drest_100'};
mol1={'POLR2A_IDR1_15_repeats'};
mol2={'MED15_IDR2'}
mytitles={'POLR2A IDR1 15 repeats - MED15 IDR2'};

Ts=[340];
reps=[1 2 3 4 5 6 7 8 9 10]; % 4, 7 have more contacts
%reps=[9 10];
%reps=[204 205 206];
frcreps=[1 2 3];

[hs, seqs]=fastaread('../POLR_idrs.fasta');


types{1}={'FWY'};
types{2}={'Q'};
types{3}={'P'};
mycolor3=[[245 130 32]/255; [31 130 65]/255; [124 40 128]/255];

colorhex={'#008181','#FEE011','#B4DCB6','#E51B4C','#F9BEBE','#25266B','#863D97','#808080','#4E64AF','#73CDD8','#BCD631','#F58230','#7F8133','#D9BDDB','#BB55A0','#9A6427'};
for i=1:length(colorhex)
    mycolor(i,:)=sscanf(colorhex{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
mycolor2=mycolor(end:-1:1,:);



%figure;
for c=1:length(constructs)
    pos=find(strcmp(hs,mol1{c})==1);
    myseq1=seqs{pos};
    pos=find(strcmp(hs,mol2{c})==1);
    myseq2=seqs{pos};

    % Get position of charged residues for seq1
    chargepos1=zeros(1,length(myseq1));
    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        for j=1:length(tmp)
            pos=strfind(myseq1,tmp(j));
            chargepos1(pos)=i;
        end
    end

    % Get position of charged residues for seq2
    chargepos2=zeros(1,length(myseq2));
    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        for j=1:length(tmp)
            pos=strfind(myseq2,tmp(j));
            chargepos2(pos)=i;
        end
    end

    count2=0;
    % Get full hamiltonian distance map
    for r=1:length(reps)
        da=load(['../' constructs{c} '/340/' num2str(reps(r)) '/inter_distanceMap_full.csv']); 
	dmapr(:,:,r)=da;
	clear da; 
    end
    mdmap=mean(dmapr,3);

    % Get FRC distance map
    for r=1:length(frcreps)
        da=load(['../' constructs{c} '/FRC_test/epsval_0.0001_epsexp_6/' num2str(frcreps(r)) '/inter_distanceMap_full.csv']); 

	dmapfrcr(:,:,r)=da;
	clear da; 
    end
    mdmapfrc=mean(dmapfrcr,3);
    
    smap=mdmap./mdmapfrc;
    
    figure;
    subplot(4,4,[6 7 8 10 11 12 14 15 16])
    matones=ones(length(myseq1),length(myseq2));
    %imagesc(smap-matones); hold on; 
    %colormap(bluewhitered)
    %caxis([-0.5 0.4])
    imagesc(smap); hold on; 
    colormap(jet)
    caxis([0.3 1.5])
    %colorbar
    title(constructs{c})

    for i=1:length(types)
        tmp=types{i};
        tmp=tmp{1};
        if i==1
            for j=1:length(tmp)
                pos=strfind(myseq1,tmp(j));
                if isempty(pos)==0
                plot(size(smap,2),pos,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
                end
            end
        elseif i==2
            for j=1:length(tmp)
                pos=strfind(myseq2,tmp(j));
		if isempty(pos)==0
                plot(pos,size(smap,1),'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
		end
            end
        elseif i==3
            for j=1:length(tmp)
                pos=strfind(myseq2,tmp(j));
		if isempty(pos)==0
                plot(pos,1,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
		end
                pos=strfind(myseq1,tmp(j));
                if isempty(pos)==0
                plot(1,pos,'o','markerfacecolor',mycolor3(i,:),'markeredgecolor','k','markersize',6); hold on;
		end
            end
        end
    end	

    %% plot net charge per residue for seq1
    bsz=5;
    for i=1:length(myseq1)-bsz+1
        ncpr(i)=(length(find(chargepos1(i:i+bsz-1)==1))-length(find(chargepos1(i:i+bsz-1)==3)))/bsz;
    end
    for i=1:length(myseq1)
        if i>=bsz & i<=length(myseq1)-bsz
            mncpr(i)=mean(ncpr(i-bsz+1:i));
        elseif i==1
            mncpr(i)=ncpr(i);
        elseif i<bsz
            mncpr(i)=mean(ncpr(1:i));
        elseif i>length(myseq1)-bsz
            tmp=length(myseq1)-i;
            mncpr(i)=mean(ncpr(end-tmp:end));
        end
    end

    subplot(4,4,[5 9 13])
    plot(1:1:length(myseq1),mncpr);
    % Extract positive and negative part
    yp = (mncpr + abs(mncpr))/2;
    yn = (mncpr - abs(mncpr))/2;
    % Plot the data using area function
    area(1:1:length(myseq1),yp,'FaceColor',mycolor3(1,:))
    hold on
    area(1:1:length(myseq1),yn,'FaceColor',mycolor3(3,:))
    xlim([1 length(myseq1)]); 
    ylim([-1 1]); 
    view(90,90)
    clear ncpr; clear mncpr; clear yp; clear yn; 

    %% plot net charge per residue for seq2
    bsz=5;
    for i=1:length(myseq2)-bsz+1
        ncpr(i)=(length(find(chargepos2(i:i+bsz-1)==2))-length(find(chargepos2(i:i+bsz-1)==3)))/bsz;
    end
    for i=1:length(myseq2)
        if i>=bsz & i<=length(myseq2)-bsz
            mncpr(i)=mean(ncpr(i-bsz+1:i));
        elseif i==1
            mncpr(i)=ncpr(i);
        elseif i<bsz
            mncpr(i)=mean(ncpr(1:i));
        elseif i>length(myseq2)-bsz
            tmp=length(myseq2)-i;
            mncpr(i)=mean(ncpr(end-tmp:end));
        end
    end

    subplot(4,4,2:4)
    plot(1:1:length(myseq2),mncpr);
    % Extract positive and negative part
    yp = (mncpr + abs(mncpr))/2;
    yn = (mncpr - abs(mncpr))/2;
    % Plot the data using area function
    area(1:1:length(myseq2),yp,'FaceColor',mycolor3(2,:))
    hold on
    area(1:1:length(myseq2),yn,'FaceColor',mycolor3(3,:))
    xlim([1 length(myseq2)]);  
    ylim([-1 1]); 

    clear smapr; clear smap; clear ncpr; clear mncpr; clear yp; clear yn; clear chargepos1; 
end

%print -painters -depsc 'two_molecule_POLR2A_IDR1_15_repeats_w_MED15_IDR2_scaled_distance_map_drest_100.eps'
%print -painters -depsc 'two_molecule_POLR2A_IDR1_15_repeats_w_MED15_IDR2_scaled_distance_map_w_colorbar_drest_100.eps'

