function save_ADCs(seq_folder, seq_name, invert_polarity, D0, R_array, N0pad, dt, restriction_geometry)

    isFEXI = contains(seq_name,{'FEXI_'});

    if invert_polarity
        seq_name = [seq_name '_np'];
    end

    load(fullfile(seq_folder, seq_name))


    %     lim = round(seq.n*1);
    %     figure(1),clf
    %     subplot(3,1,1)
    %     plot(seq.gc,'o-')
    %     xlim([0 lim])
    %     subplot(3,1,2)
    %     plot(seq.b,'o-')
    %     xlim([0 lim])
    %     subplot(3,1,3)
    %     plot(seq.tm,'o-')
    %     xlim([0 lim])
    %     return

    display(sprintf('%s: delta max = %g ms', seq_name, seq.delta_max*1e3));

    delta_max = seq.delta_max; % determines transverse time (must be >= max delta )

    Tmax = 2*(2*delta_max + seq.delta_pi) + 2*seq.delta_c + max(seq.tm_array);

    Nmax = round(Tmax/dt);
    NFT = Nmax + N0pad;
    f = linspace(-1,1,NFT)'/2/dt;

    flim = 1000*[0 1];
    ind_PS = find(f >= flim(1) & f<= flim(2));
    fPS = f(ind_PS);

    switch restriction_geometry
        case 1
            Dw = fDwSphere((2*pi*fPS).^2,R_array,D0,0,50)'; % spherical restriction
        case 2
            Dw = fDwCylinder((2*pi*fPS).^2,R_array,D0,0,50)'; % cylindrical restriction
        case 3
            Dw = fDwPlane((2*pi*fPS).^2,R_array,D0,0,50)'; % planar restriction
    end


    Seq = seq; % required by parfor?

    tic
    ADC = zeros(Seq.n,length(R_array));
    %     for c = 1:Seq.n
    parfor c = 1:Seq.n
        if isFEXI
            [~, q] = def_FEXI(Seq.delta, Seq.delta_max, Seq.g1(c), Seq.g2(c), Seq.tm(c), Seq.gc(c), Seq.delta_c, Seq.delta_pi, dt);
        else
            %[~, q] = def_sequence(Seq.delta(c), delta_max, Seq.g(c), Seq.tm(c), Seq.gc(c), Seq.delta_c, Seq.delta_pi, dt);
            [~, q] = def_sequence(Seq.delta(c), delta_max, Seq.g1(c), Seq.g2(c), Seq.tm(c), Seq.gc(c), Seq.delta_c, Seq.delta_pi, dt);

        end
        %                     figure(1),clf
        %                     subplot(2,1,1), plot(g)
        %                     subplot(2,1,2), plot(q)

        s = power_spectrum(q, dt, NFT, ind_PS);
        ADC(c,:) = Dw * s;

        %             figure(1),clf
        %             plot(s(1:end/8))
    end
    toc

    DvsR.D0 = D0;
    DvsR.ADC = ADC;
    DvsR.R_array = R_array;

    path = fullfile(seq_folder, [seq_name '_DvsR']);
    save(path,'DvsR')

