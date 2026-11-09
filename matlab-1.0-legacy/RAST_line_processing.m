function [sol]=RAST_line_processing(Euler,I,ii,no,path,coefmat,degrot,error,cuth,cutl,nei_deg,pixsize,numlines,scorefigs,crystal)

% Original version by Vincent Gagneur
%
%This function processes every crystal orientation (Euler), picture (I),
%combinations. It finds what each detected line matches in terms of
%theoretical slip plane type slip traces. This code also gives various
%statistics (match counts), generates and tests randomly oriented lines, as
%well as generates figures.
%
% By default only two "special conditions" are checked for for figure
% making, lines fitting these conditions show up in magenta: 
% - in bcc if a line is unambiguously of <111> type i.e. it only
% matche 112 and or 123 orientation
% - in fcc if a line is ambiguously oriented (as only 4 slip planes, this
% is significantly less frequent than in bcc)
%
% NOTE: in this code nei = score threshold ST in the paper

[h,l]=size(I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gets Euler and calculates theoretical slips with Euler_to_Slip

[slip100,slip110,slip112,slip123]=RAST_euler_to_slip(Euler,crystal); % <-- Euler_to_Slip.m must be in the same folder as Slip_Trace_Analyser_single.m Calculates theoretical slip traces for each bcc plane

%normalisation of vectors
slips = {slip100, slip110, slip112, slip123};

for j=1:length(slips)
    slipette=slips{j}; 
    for i=1:length(slipette)
        slipette(i,:)=slipette(i,:)/norm(slipette(i,:));
    end
    slips{j}=slipette;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The next loop converts the vectors (slip orientations) to corresponding slopes

CO={}; %matrix where corresponding slopes will be stored
for i = 1:length(slips)
    slip=transpose(slips{i});
    coefstore=zeros(1,length(slip));
    for j = 1:length(slip)
        slope=(slip(2,j))/(-slip(1,j)); %determines the slope of each slip, slope = y / x, correction x --> -x warning can change if euler does not come  from  mtex
        coefstore(j)=slope;
    end
    CO{i}=coefstore;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Slip analysis loop 

%g = waitbar(0,sprintf('Processing every line detected (picture %.0f of %.0f)',[ii no]));
fprintf('Processing every line detected (picture %.0f of %.0f)',[ii no])
fprintf('\n');


count=0;

best_all=zeros(cuth,cutl,numlines,4); %will store best match
second_all=zeros(cuth,cutl,numlines,4); % will store secondary matches
matches_all=zeros(cuth,cutl,numlines,4); % will store both

for i = 1:1:(cuth)%
    for j = 1:1:(cutl)

    A_1=coefmat(i,j,numlines+1);

    %next loop compares each detected line to the theoretical plane
    %orientations, and saves what plane types are matched within IR (error) degrees
    if A_1~=0 %not on a masked area

        for ln=1:1:numlines

            co=coefmat(i,j,ln);

            count=count+1;
            %waitbar(count/(cuth*cutl))

            [best,second,matches]=RAST_slip_trace_analyser_single(CO,slips,co,error); % calls slip matching code
    
            %store best matches info for figure making
            if best(1,1,1,1)==1
                bestplane=100; %stores as double because easier
            elseif best(1,1,1,2)==1
                bestplane=110;
            elseif best(1,1,1,3)==1
                bestplane=112;
            elseif best(1,1,1,4)==1
                bestplane=123;
            else 
                bestplane=0;
            end

            coefmat(i,j,numlines+6+6*(ln-1))=bestplane; %for figure 

            best_all(i,j,ln,:)=best; 
            second_all(i,j,ln,:)=second;
            matches_all(i,j,ln,:)=matches;

        end

    end

    end
end

%close(g)

%Now the next loop will count the total amount of each match, unique match,
%double or triple match. Not counting quadruple match because we learn nothing from their analysis? (not sure why I did not code it).
for nei=0:1:(20) %nei is the score threshold S, here goes from 0 to 20

    m100=0;
    m110=0;
    m112=0;
    m123=0;

    u100=0;
    u110=0;
    u112=0;
    u123=0;

    u111=0;

    m100110=0;
    m100112=0;
    m100123=0;
    m110112=0;
    m110123=0;
    m112123=0;

    m110112123=0;
    m100112123=0;
    m100110123=0;
    m100110112=0;

    nomatch=0;

    lcount=0;

    for i = 1:1:(cuth)
        for j = 1:1:(cutl)
            A_1=coefmat(i,j,numlines+1);

            if A_1~=0 %not on a masked area

                for ln=1:1:numlines %numlines lines per segment
            
                    if coefmat(i,j,numlines+5+6*(ln-1))>=nei % checks if line score is >= to score threshold ST (nei)
            
                        lcount=lcount+1; %simply counting total number of lines
                        
                        %recover all matches:
                       
                        match=matches_all(i,j,ln,:);

                        if match(1,1,1,1)==1
                            m100=m100+1; %counts how many times each plane was matched
                        end

                        if match(1,1,1,2)==1
                            m110=m110+1;
                        end

                        if match(1,1,1,3)==1
                            m112=m112+1;
                        end

                        if match(1,1,1,4)==1
                            m123=m123+1;
                        end
                        
                        %counting if matching nothing

                        if  sum(match(1,1,1,:))==0
                            nomatch=nomatch+1;
                        end

                        %now counting unambiguously <111> slip for fun if
                        %bcc

                        if crystal==1
                            if  sum(match(1,1,1,:))>0
                                if (match(1,1,1,1)==1)||(match(1,1,1,2)==1) 
                                    %nothing
                                else % = it has matched something, but not 100 or 110
                                    u111=u111+1;
                                    coefmat(i,j,numlines+6+6*(ln-1))=111; %for figure making, line will be magenta
                                end
                            end
                        end

                        % Now couting how many unique matches

                        if sum(match(1,1,1,:))==1

                            if match(1,1,1,1)==1
                                u100=u100+1;
                            end

                            if match(1,1,1,2)==1
                                u110=u110+1;
                            end

                            if match(1,1,1,3)==1
                                u112=u112+1;
                            end

                            if match(1,1,1,4)==1
                                u123=u123+1;
                            end

                        % if analysing fcc, may be interesting to know on figure
                        % which lines are ambiguously oriented as it is much
                        % less frequent than in bcc. Can also be deduced
                        % from the stats (simply count total amount of
                        % unique match/total amount of lines):

                        elseif crystal==2
                            if sum(match(1,1,1,:))>0 %means neither 0 nor 1
                            coefmat(i,j,numlines+6+6*(ln-1))=111; %for figure making, line will be magenta, matche something but also something else.
                            end
                        end

                        %now counting each double match

                        if sum(match(1,1,1,:))==2

                            if match(1,1,1,1)==1 && match(1,1,1,2)==1
                               m100110=m100110+1;
                            end
                            if match(1,1,1,1)==1 && match(1,1,1,3)==1
                               m100112=m100112+1;
                            end
                            if match(1,1,1,1)==1 && match(1,1,1,4)==1
                               m100123=m100123+1;
                            end

                            if match(1,1,1,2)==1 && match(1,1,1,3)==1
                               m110112=m110112+1;
                            end
                            if match(1,1,1,2)==1 && match(1,1,1,4)==1
                               m110123=m110123+1;
                            end

                            if match(1,1,1,3)==1 && match(1,1,1,4)==1
                               m112123=m112123+1;
                            end

                        end

                        %now counting each triple match

                        if sum(match(1,1,1,:))==3

                            if match(1,1,1,1)==0 
                               m110112123=m110112123+1;
                            end
                            if match(1,1,1,2)==0 
                               m100112123=m100112123+1;
                            end
                            if match(1,1,1,3)==0 
                               m100110123=m100110123+1;
                            end
                            if match(1,1,1,4)==0 
                               m100110112=m100110112+1;
                            end

                        end

                    end    
                           
                 end
                      
           end
        end
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Slip analysis loop for randomly generated lines, also called "test" lines (basically same loop)

    tm100=0;
    tm110=0;
    tm112=0;
    tm123=0;

    tu100=0;
    tu110=0;
    tu112=0;
    tu123=0;

    tu111=0;

    tm100110=0;
    tm100112=0;
    tm100123=0;
    tm110112=0;
    tm110123=0;
    tm112123=0;

    tm110112123=0;
    tm100112123=0;
    tm100110123=0;
    tm100110112=0;

    tnomatch=0;

    tcount=0;

    for i = 0.25:0.25:180 % lines of orientations 0.25 to 180° are generated by the code, to evaluate random line orientations likelihood of matching each condition.

        tcount=tcount+1;

        [best,second,match]=RAST_slip_trace_analyser_single(CO,slips,i,error); % calls slip matching code

        if match(1,1,1,1)==1
            tm100=tm100+1; %counts how many times each plane was matched
        end

        if match(1,1,1,2)==1
            tm110=tm110+1;
        end

        if match(1,1,1,3)==1
            tm112=tm112+1;
        end

        if match(1,1,1,4)==1
            tm123=tm123+1;
        end
        
        %counting if matching nothing

        if  sum(match(1,1,1,:))==0
            tnomatch=tnomatch+1;
        end

        %now counting unambiguously <111> slip for fun

        if  sum(match(1,1,1,:))>0
            if (match(1,1,1,1)==1)||(match(1,1,1,2)==1)
                %nothing
            else
                tu111=tu111+1;
            end
        end

        % Now counting how many unique matches

        if sum(match(1,1,1,:))==1

            if match(1,1,1,1)==1
                tu100=tu100+1;
            end

            if match(1,1,1,2)==1
                tu110=tu110+1;
            end

            if match(1,1,1,3)==1
                tu112=tu112+1;
            end

            if match(1,1,1,4)==1
                tu123=tu123+1;
            end

        end

        %now counting each double match

        if sum(match(1,1,1,:))==2

            if match(1,1,1,1)==1 && match(1,1,1,2)==1
               tm100110=tm100110+1;
            end
            if match(1,1,1,1)==1 && match(1,1,1,3)==1
               tm100112=tm100112+1;
            end
            if match(1,1,1,1)==1 && match(1,1,1,4)==1
               tm100123=tm100123+1;
            end

            if match(1,1,1,2)==1 && match(1,1,1,3)==1
               tm110112=tm110112+1;
            end
            if match(1,1,1,2)==1 && match(1,1,1,4)==1
               tm110123=tm110123+1;
            end

            if match(1,1,1,3)==1 && match(1,1,1,4)==1
               tm112123=tm112123+1;
            end

        end

        %now counting each triple match

        if sum(match(1,1,1,:))==3

            if match(1,1,1,1)==0 
               tm110112123=tm110112123+1;
            end
            if match(1,1,1,2)==0 
               tm100112123=tm100112123+1;
            end
            if match(1,1,1,3)==0 
               tm100110123=tm100110123+1;
            end
            if match(1,1,1,4)==0 
               tm100110112=tm100110112+1;
            end

        end
    
    end



    nlines=lcount; %calcualting total number of lines detected
    tnlines=tcount; %same but for test lines generated

    %now calculating corresponding percentages "p" for each value.

    pm100=m100/nlines*100;
    pm110=m110/nlines*100;
    pm112=m112/nlines*100;
    pm123=m123/nlines*100;
    pu100=u100/nlines*100;
    pu110=u110/nlines*100;
    pu112=u112/nlines*100;
    pu123=u123/nlines*100;
    pu111=u111/nlines*100;
    pm100110=m100110/nlines*100;
    pm100112=m100112/nlines*100;
    pm100123=m100123/nlines*100;
    pm110112=m110112/nlines*100;
    pm110123=m110123/nlines*100;
    pm112123=m112123/nlines*100;
    pm110112123=m110112123/nlines*100;
    pm100112123=m100112123/nlines*100;
    pm100110123=m100110123/nlines*100;
    pm100110112=m100110112/nlines*100;

    pnomatch=nomatch/nlines*100;

    %now calculating corresponding percentages "p" for each value, for test lines "t".

    ptm100=tm100/tcount*100;
    ptm110=tm110/tcount*100;
    ptm112=tm112/tcount*100;
    ptm123=tm123/tcount*100;
    ptu100=tu100/tcount*100;
    ptu110=tu110/tcount*100;
    ptu112=tu112/tcount*100;
    ptu123=tu123/tcount*100;
    ptu111=tu111/tcount*100;
    ptm100110=tm100110/tcount*100;
    ptm100112=tm100112/tcount*100;
    ptm100123=tm100123/tcount*100;
    ptm110112=tm110112/tcount*100;
    ptm110123=tm110123/tcount*100;
    ptm112123=tm112123/tcount*100;
    ptm110112123=tm110112123/tcount*100;
    ptm100112123=tm100112123/tcount*100;
    ptm100110123=tm100110123/tcount*100;
    ptm100110112=tm100110112/tcount*100;

    ptnomatch=tnomatch/tcount*100;


    %and returning all values for current ST (nei):

    y=[m100 m110 m112 m123 u100 u110 u112 u123 u111...
    m100110 m100112 m100123 m110112 m110123 m112123...
    m110112123 m100112123 m100110123 m100110112 nomatch...
    tm100 tm110 tm112 tm123 tu100 tu110 tu112 tu123 tu111...
    tm100110 tm100112 tm100123 tm110112 tm110123 tm112123...
    tm110112123 tm100112123 tm100110123 tm100110112 tnomatch...
    pm100 pm110 pm112 pm123 pu100 pu110 pu112 pu123 pu111...
    pm100110 pm100112 pm100123 pm110112 pm110123 pm112123...
    pm110112123 pm100112123 pm100110123 pm100110112 pnomatch...
    ptm100 ptm110 ptm112 ptm123 ptu100 ptu110 ptu112 ptu123 ptu111...
    ptm100110 ptm100112 ptm100123 ptm110112 ptm110123 ptm112123...
    ptm110112123 ptm100112123 ptm100110123 ptm100110112 ptnomatch...
    nlines tnlines];
    
    sol{nei+1}=y;

    %%%%%%%%%%%%%%%%%
    %%%
    %figure plotting the lines selected (current nei condition)

    bnlines=sum(y(1:5));

    if bnlines~=0 && ismember(nei,scorefigs)==1 % <--- Plotting figures can be slow, so you can choose which score ST to plot figures for. Note that some lines can have negative scores, and therefore ST=0 already deletes parts of the detected lines.

        RAST_Figuremaking(I,ii,nei,pixsize,error,cuth,cutl,numlines,coefmat,Euler,degrot,path,crystal)

    end

    end
end

