% Converts 20x20 pixel matrix to PBS structure

function er = pixelToEr(pbs, er, erSi, Ny, Nx, dy, dx, pbs_size)
    Npix = 20;
    pix_size = 120e-9;
%     er(66:187,66:187) = 12;
%     er(151:173,186:end) = 12;
%     er(79:101,186:end) = 12;
%     er(115:137,:) = 12;
    origin = [(Ny/2)+1 (Nx/2)+1];
    btmleftcornerY = origin(1)-pbs_size/2/dy;
    btmleftcornerX = origin(2)-pbs_size/2/dx;
    pixEr = pix_size/dx;
    for i=1:Npix
        p1 = btmleftcornerY;
        for j=1:Npix
            p2 = floor(btmleftcornerY+(j)*pixEr)-1;
            if pbs(i,j) == 1
                er(p1:p2,...
                   floor(btmleftcornerX+(i-1)*pixEr):floor(btmleftcornerX+(i)*pixEr)-1) = erSi;
            end
            p1=p2+1;
        end
    end
end

