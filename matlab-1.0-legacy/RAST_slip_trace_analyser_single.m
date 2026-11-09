function  [best,second,matches]=RAST_slip_trace_analyser_single(CO,slips,lineslope,error)

%%%%%%%%%%%%%%
% Automatic matching of theoretical slip traces vs detected lines, by Vincent Gagneur
% 
% The use of this code requires another code "Euler_to_Slip.m" which
% calculates the theoretical slip traces based the Euler angles of a crystal orientation.
%
% Currently part of the automated loop, need a line orientation and eulers
% angles as entries, as well as the identification range "error" (IR in paper) which
% is the maximal angular difference between theoretical slip orientations and
% detected line orientations to consider plane type as matching the line (i.e. as being possibly the originating plane type of the observed slip).
% best and secondary matches come out. One line per loop.
%
%
%
% Note that the code also calculates the angular difference between detected
% lines and theoretical orientation, which is only used for finding the best match for figure making, but unusued for the statistical
% analysis. Further code expansion could easily give access to these values if of interest.

%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Gets Euler and calculates theoretical slips with Euler_to_Slip
% 

second=zeros(1,1,1,4);
best=zeros(1,1,1,4);

% 
% [slip100,slip110,slip112,slip123]=RAST_euler_to_slip(Euler,crystal); % <-- Euler_to_Slip.m must be in the same folder as Slip_Trace_Analyser_single.m Calculates theoretical slip traces for each bcc plane
% 
% %normalisation of vectors
% slips = {slip100, slip110, slip112, slip123};
% 
% for j=1:length(slips)
%     slipette=slips{j}; 
%     for i=1:length(slipette)
%         slipette(i,:)=slipette(i,:)/norm(slipette(i,:));
%     end
%     slips{j}=slipette;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %The next loop converts the vectors (slip orientations) to corresponding slopes
% 
% CO={}; %matrix where corresponding slopes will be stored
% for i = 1:length(slips)
%     slip=transpose(slips{i});
%     coefstore=zeros(1,length(slip));
%     for j = 1:length(slip)
%         slope=(slip(2,j))/(-slip(1,j)); %determines the slope of each slip, slope = y / x, correction x --> -x warning can change if euler does not come  from  mtex
%         coefstore(j)=slope;
%     end
%     CO{i}=coefstore;
% end

lineslope=tand(lineslope); %need converting to make it same format as theoretical slips

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparing the slopes of expected slip traces and the line analysed


slip100co=[];
slip110co=[];
slip112co=[];
slip123co=[];
degoff100=[];
degoff110=[];
degoff112=[];
degoff123=[];


%Main comparison loop:
for i = 1:length(CO)
    slipsys=CO{i};
    matched=0;
    sloped=zeros(1);
    degoff=zeros(1);
    for j = 1:length(slipsys)
        slope=slipsys(j);
        refdeg=atand(slope);
        linedeg=atand(lineslope);

        if (((refdeg-error)<(linedeg)) & ((linedeg)<(refdeg+error)))||(((refdeg-error+180)<(linedeg)) & ((linedeg)<(refdeg+error+180)))||(((refdeg-error-180)<(linedeg)) & ((linedeg)<(refdeg+error-180)))
            matched=matched+1; 
            
            sloped(numel(zeros(1,matched)))=0;
            sloped(1,matched)=slope; %stores line as matching line
            
            deg=abs(atand(slope)-atand(lineslope));
            
            degoff(numel(zeros(1,matched)))=0;
            degoff(1,matched)=deg; % also stores the error of the matched line NB: currently unused
        end
    end

    %for debugging + stores the slip planes that matched for each
    %slip system:

    if matched~=0
     if i==1
            %disp(strcat('slip100 was matched ',int2str(matched),'times'))
            slip100co=sloped;
            degoff100=degoff;
     elseif i==2
            %disp(strcat('slip110 was matched ',int2str(matched),'times'))
            slip110co=sloped;
            degoff110=degoff;
     elseif i==3
            %disp(strcat('slip112 was matched ',int2str(matched),'times'))
            slip112co=sloped;
            degoff112=degoff;
     elseif i==4
            %disp(strcat('slip123 was matched ',int2str(matched),'times'))
            slip123co=sloped;
            degoff123=degoff;
     end
   end

end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % match storing procedure, also finds which match is the most closely
    % oriented
    
    degoffmat={degoff100,degoff110,degoff112,degoff123};
    
    co=[slip100co,slip110co,slip112co,slip123co];
    de=[degoff100,degoff110,degoff112,degoff123];
    
    for i = 1:length(slips)
        slip=transpose(slips{i});

        slipsys=CO{i};
        degoffsys=degoffmat{i};
        Y=slipsys;

        hasmatched=0;
        bestmatch=0;
        matchissime=0;
        
        beststring='NO MATCH';
        names={'(100)','(110)','(112)','(123)'};
        
        
        
        for j = 1:length(slip) 
            if ismember(Y(j),co)==1    
               hasmatched=1;
            else

            end
            if hasmatched==1
                for z = 1:length(co)
                    if Y(j)==co(z)
                        x=z;
                    end
                end
                if de(x)==min(degoffsys) 
                   bestmatch=1; %nb variable is called bestmatch, because it is the best match amongst the different possible planes within one specific plane type =/= best matching plane type
                   if de(x)==min(de) 
                       matchissime=1; % if smallest angular difference of all plane types, it is the best matching plane type
                   end
                end
            end
        end
        %Stored infos of which match are best or secondary
        if matchissime==1
             best(1,1,1,i)=1;
        elseif bestmatch==1
             second(1,1,1,i)=1;
        end

    end
    matches=best+second;
end
            
            
            
            
            
            


            
    