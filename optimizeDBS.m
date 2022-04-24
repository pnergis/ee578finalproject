% Direct Binary Search Optimizer

function [FOM_TE, FOM_TM, pbs] = optimizeDBS(pbs)
   
    % initialize FOM for TE and TM
    FOM_TE_old = 0;
    FOM_TM_old = 0;
    % stores pixels changed (400 max pixels b/c 20x20)
    history = zeros(400,2);
    % count # of iterations
    i=1;
    % preallocate for 1000 iterations, keep track of improvement over
    % iterations
    FOM_TE = zeros(1000,1);
    FOM_TM = zeros(1000,1);
    % begin optimization loop which ends once TE and TM have minimum 75%
    % transmission efficiency each in their respective output waveguides
    while (FOM_TE_old<85 || FOM_TM_old<85)
        % run FDTD for new FOMs
        FOM_TM(i) = FDTD_2D_TM(pbs,0,0,0);
        FOM_TE(i) = FDTD_2D_TE(pbs,0,0,0);
        
        % if results worsened for either TE and TM
        if FOM_TE(i)<FOM_TE_old || FOM_TM(i)<FOM_TM_old
                % unflip the pixel
            if pbs(nxtpix(1),nxtpix(2)) == 1
                pbs(nxtpix(1),nxtpix(2)) = 0;
            else
                pbs(nxtpix(1),nxtpix(2)) = 1;
            end
        else
            % if pixel flip improved result, show new design
            figure (1)
            pcolor(pbs)
            shading interp
            xlabel('x')
            ylabel('y') 
            title('\epsilon TM');
            set(gca,'YDir','normal')
        end
        disp([i FOM_TM(i) FOM_TE(i)])
        i=i+1; % update iteration count
        %check if other flips are possible
        if nnz(history) ~= 800 % 2*400 possible pixels * 2 because (m,n) storage
            FOM_TE_old = FOM_TE;
            FOM_TM_old = FOM_TM;
            nxtpix = [randi([1 20],1) randi([1 20],1)];
            checkNew = ismember(history(:,1:2),nxtpix,'rows');
            % keep randomly changing until nxtpix is a new pixel to change
            while nnz(checkNew)~=0
                nxtpix = [randi([1 20],1) randi([1 20],1)];
                checkNew = ismember(history(:,1:2),nxtpix,'rows');
            end
            % update pbs pixels
            if pbs(nxtpix(1),nxtpix(2)) == 1
                pbs(nxtpix(1),nxtpix(2)) = 0;
            else
                pbs(nxtpix(1),nxtpix(2)) = 1;
            end
            % update history of changed pixels
            history(i,:) = nxtpix;
        else
            % clear the history list and rerun through pixels with new
            % optimized initial pbs until desired FOM achieved
            history = zeros(400,2);
        end
    end
end