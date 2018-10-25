function [estMFPT] = RunOU2D_MTarg_WtCorr(c,MaxIter,MFPTIter)
%Model Params--User input
tauSlow=8.18E-5;
tauFast=5.22E-7;
sigmax=4.27;
sigmay=4.27;
dt=0.0005*tauFast;
%c=sqrt(0.480);


%WE Simulation Params--User input
tau_WE=1E-8;   %WE iteration length
MTarg=200;      %Target number of replicas per bin
%MaxIter=10000;  %maximum number of tau_WE iterations
PlotIter=MaxIter; %number of iterations to keep for plotting results (takes from the end)
%MFPTIter=9000;  %number of iterations to keep for computing MFPT
InitialConfigFile=1; %1 to read initial replicas from file, 0 to generate initial config
StoreSz=3;
RepArSz=2; %size of replica array, i.e., the size of state configuration
WtTot=1E4; %initial value of the total weight
SaveEvery=5E2;
counter=0;
fluxfilename=['flux_' num2str(c)];

%Setting various parameters
taufold=round(tau_WE/dt); %this is number of dt-steps in each iter. tau_WE is taufold*dt
TargetX=25;
BinDefs=[-15,-10,-5,0,5:1:9,10:.25:20,20.25:.1:TargetX]; %Bin Assignments in terms of 1D Order Param
NBins=numel(BinDefs)+1;
NReps=NBins*MTarg; %Total number of replicas in simulation
keep_bin_weights=zeros(NBins,MaxIter);

fluxes=zeros(MaxIter,1);    %initialize array of fluxes
if InitialConfigFile %loading stored replicas/weights from file
    load replicas
    load weights
    weights=weights/sum(weights)*WtTot; %just in case weights don't sum to WtTot, make them!
else %Generate a random initial configuration for simulation
    initialpos=randn(2,1)*sigmax;
    %make sure the initial replica is not in target state, if so, regenerate
    while numel(ReachedTarget(initialpos))>0
        initialpos=randn(2,1)*sigmax;
    end
    replicas=initialpos
    weights=1*WtTot;
end


for WEiter=1:MaxIter  %loop over tau-steps
    WEiter
    
    if (WEiter-counter*SaveEvery)/SaveEvery >= 1
        save WEiter WEiter
        counter=counter+1;
        save(fluxfilename,'fluxes')
    end
    fluxes(WEiter)=0; %currently not initializing the fluxes array
    %propagate Brownian Dynamics
    for BDiter=1:taufold
        
        replicas=PropagateDynamics(replicas);
        
        %         %implement "stop" condition and reintroduction: check if
        %         evacuation has occurred and transfer weight from evacuated
        %         replicas
        
        evac_inds=ReachedTarget(replicas);
        if numel(evac_inds)>0
            transfer_wt=weights(evac_inds); %weight of evacuated replicas--to be reintroduced
            weights(evac_inds)=0; %evacuated replicas have weight set to 0
            fluxes(WEiter)=fluxes(WEiter)+sum(transfer_wt)/WtTot; %add evacuated weights to fluxes at each iter.
            %remove the evacuated replica from the system
            replicas=replicas(:,weights>0);
            weights=weights(weights>0);
            %the evacuated weights are reintroduced by renormalizing remaining weight to 1
            weights=weights/sum(weights)*WtTot;
        else
            transfer_wt=0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %WE step--return update_weights and update_replicas from weights and
    %replicas. (Also bin_weights: the summed weights in each bin)
    
    %weight correction - to correct any potential loss/addition of weight
    weights=weights/sum(weights);
    
    [update_replicas, update_weights, bin_weights] = WE_Step_MTarg(replicas,weights);
    
    keep_bin_weights(:,WEiter)=bin_weights;
    
    %remove the padding 0's from replicas and weights to continue
    %propagating dynamics
    getinds=find(update_weights>0);
    replicas=update_replicas(:,getinds);
    weights=update_weights(update_weights>0);
end
stind=max([1,WEiter-MFPTIter]);
estFlux=mean(fluxes(stind:WEiter)/tau_WE);
estMFPT=1/estFlux
PlotResults(keep_bin_weights,fluxes)

save keep_bin_weights keep_bin_weights
save fluxes fluxes
save replicas replicas
save weights weights

    function newreplicas=PropagateDynamics(replicas)
        %OU process
        newreplicas=replicas - 1./[tauSlow; tauFast].*replicas*dt + ([sigmax; sigmay]./(sqrt([tauSlow; tauFast])))*sqrt(2*dt).*randn(size(replicas));
    end

    function BinIndex=AssignBinIndex(replicas) %assigns a particular order param value to bin
        %Bin index should start at 1 and not exceed N
        OrderParam=Compute_OrderParam(replicas);
        OrderParam=OrderParam(:);
        CompArray=OrderParam>=BinDefs;
        BinIndex=sum(CompArray,2)+1;
    end

    function Indices=ReachedTarget(replicas) %finds out if any replicas are in target state
        OrderParam=Compute_OrderParam(replicas);
        Indices=find(OrderParam>TargetX);
    end

    function OrderParam=Compute_OrderParam(replicas) %compute Order Param for each replica
        OrderParam=c*replicas(1,:)+sqrt(1-c^2)*replicas(2,:);
    end


    function [update_replicas, update_weights, bin_weights]=WE_Step_MTarg(replicas,weights)
        %Determine which bins the replicas are in after tau_WE
        BinIndex=AssignBinIndex(replicas); %list of bin assignments for each replica
        bin_weights=zeros(NBins,1);
        update_replicas=zeros([RepArSz],MTarg*NBins);
        update_weights=zeros(MTarg*NBins,1);
        
        for binloop=1:NBins %loop over each bin for the splitting/merging step
            
            update_replicas_bin=zeros(RepArSz,MTarg);  %initialize empty array of updated replicas to populate in WE step
            update_weights_bin=zeros(MTarg,2);  %initialize empty array of updated weights to populate in WE step
            fixed_bininds=(binloop-1)*(MTarg)+[1:MTarg]; %the fixed indices where the replicas will be placed after updating
            
            repinds_bin=find(BinIndex==binloop); %find indices of replicas in current bin
            replicas_currbin=replicas(:,repinds_bin); %just the replicas currently in bin
            weights_currbin=weights(repinds_bin); %weights in current bin
            bin_weights(binloop)=sum(weights_currbin); %summed weight in current bin
            
            M=nnz(weights_currbin); %number of replicas currently in this bin
            
            if M>0
                update_weights_bin(1:M,1)=weights_currbin;
                update_weights_bin(1:M,2)=[1:M]';
                while M<MTarg
                    %to increase number of replicas, select biggest replica and split it
                    update_weights_bin=sortrows(update_weights_bin,1,'descend');
                    halfweight=max(update_weights_bin(1,1)/2,realmin);
                    update_weights_bin(M+1,1)=halfweight;
                    update_weights_bin(M+1,2)=update_weights_bin(1,2);
                    update_weights_bin(1,1)=halfweight;
                    M=nnz(update_weights_bin(:,1));
                end
                while M>MTarg
                    %to decrease the number of replicas, select smallest pairs and merge them
                    update_weights_bin=sortrows(update_weights_bin,1,'descend');
                    smallest2=update_weights_bin(end-1:end,1);
                    %select one of the smallest 2 to keep configuration
                    lumped_wt=sum(smallest2);
                    norm_wt1=smallest2(1)./lumped_wt;
                    if rand<norm_wt1
                        update_weights_bin(end,1)=0;
                        update_weights_bin(end-1,1)=lumped_wt;
                    else
                        update_weights_bin(end-1,1)=0;
                        update_weights_bin(end,1)=lumped_wt;
                    end
                    update_weights_bin=update_weights_bin(update_weights_bin(:,1)>0,:);
                    M=nnz(update_weights_bin(:,1));
                end
                update_replicas_bin(:,1:MTarg)=replicas_currbin(:,update_weights_bin(:,2));
            else
            end
            %place the updated replcias and weights in the full array
            update_weights(fixed_bininds)=update_weights_bin(:,1); %put weights in full array
            update_replicas(:,fixed_bininds)=update_replicas_bin;
        end
    end

    function []=PlotResults(keep_bin_weights,fluxes)
        close all
        
        figure(1)
        subplot(2,1,1)
        stind=max([1,WEiter-PlotIter]);
        x=[stind:WEiter];
        semilogy(x,fluxes(stind:WEiter)/tau_WE,'+')
        y=repmat(estFlux,size(x));
        hold on
        semilogy(x,y)
        ylabel('flux')
        xlabel('iteration')
        txt=['MFPT = ' num2str(estMFPT)];
        text(stind+10,estFlux+.5*estFlux,txt)
        
        subplot(2,1,2)
        X=[BinDefs,15];
        plot(X,keep_bin_weights(:,x),'.','Color',[0.6 0.6 0.6])
        hold on
        plot(X,mean(keep_bin_weights(:,x),2),'-+r')
        xlabel('X')
        ylabel('Probability')
        FigName=['Results' num2str(c*10) '.png' ];
        %print -dpng ResultsFig
        print(FigName,'-dpng')
    end

end

