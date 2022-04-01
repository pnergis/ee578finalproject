% PBS structure

function ep = pixelToEp(pbs, Ny, Nx, dy, dx, ep, er, pbs_size)
    Npix = 20;
    pix_size = 240e-9;
    for i=1:Npix
        for j=1:Npix
            if pbs(i,j) == 1
                ep((Ny+pbs_size/dy-j*pix_size/dy):(Ny+pbs_size/dy-(j-1)*pix_size/dy),...
                    (Nx-pbs_size/dx+(i-1)*pix_size/dx):(Nx-pbs_size/dx+i*pix_size/dx)) = er;
            end
        end
    end
end
