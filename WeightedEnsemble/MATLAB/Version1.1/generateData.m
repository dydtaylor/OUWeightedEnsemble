for j = 1:10
    RunOU1D_AlwaysMTarg()
    filename = "fluxes" + j;
    save(filename, 'fluxes')
end