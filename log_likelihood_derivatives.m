function logL = log_likelihood_derivatives(Bz_list, phi_list, bz, Bperp, alpha, beta, zeta, S_exp, sigma)

    % Evaluates a Gaussian log-likelihood for PL data given a simplified NV
    % crossing-lineshape model. The model includes first-order uncertainty
    % propagation from the swept variables (Bz, phi) via finite-difference
    % derivatives of the PL signal.
    
    % Inputs
    % Bz_list : vector [nBz x 1] of longitudinal sweep values (Tesla)
    % phi_list : vector [1 x nPhi] of in-plane angles (radians)
    % bz : scalar axial field component in the laboratory frame (Tesla)
    % Bperp : scalar transverse field magnitude in the laboratory frame (Tesla)
    % alpha : scalar orientation angle (radians)
    % beta : scalar orientation angle (radians)
    % zeta : scalar orientation angle (radians)
    % S_exp : matrix [nBz x nPhi] or [nPhi x nBz] of experimental PL (unitless)
    % sigma : scalar baseline noise s.d. of PL (same units as S_exp)
    
    % Output
    % logL : scalar Gaussian log-likelihood value
    

    % Width and contrast are phenomenological and tuned to the experiment
    width = 1.5e-4;
    contrast = 0.0042;
    
    % Derivative steps (dBz, dphi) default to 1/50 of the sampling step to
    % stabilize finite-difference estimates; uncertainty in Bz and phi is
    % set to 1/5 of the sampling step.
    dBz = (Bz_list(2)-Bz_list(1))/50;
    sigma_Bz = (Bz_list(2)-Bz_list(1))/5;

    if length(phi_list)>1
        dphi = (phi_list(2)-phi_list(1))/50;
        sigma_phi = (phi_list(2)-phi_list(1))/5;
    else
        dphi = 5*pi/180/50;
        sigma_phi = 1*pi/180;
    end

    % Magnetic field components in the diamond frame
    Bx_NV = @(Bz,phi) 1/6*(sqrt(2)*(bz+Bz)*sin(beta)*(2 + sqrt(3)*sin(zeta)) ...
        + sqrt(2)*Bperp*sin(beta)*(-1 + sqrt(3)*sin(zeta))*(cos(alpha - phi) ...
        - sin(alpha - phi)) - sqrt(3)*Bperp*sin(zeta)*(cos(alpha - phi) + ...
        sin(alpha - phi)) + cos(zeta)*(sqrt(2)*(bz+Bz)*sin(beta) + ...
        Bperp*cos(alpha - phi)*(3 + sqrt(2)*sin(beta)) + ...
        Bperp*(3 - sqrt(2)*sin(beta))*sin(alpha - phi)) + ...
        cos(beta)*(-2*(bz+Bz)*(-1 + cos(zeta) + sqrt(3)*sin(zeta)) ...
        +  Bperp*cos(alpha - phi)*(2 + cos(zeta) + sqrt(3)*sin(zeta)) ...
        -   Bperp*(2 + cos(zeta) + sqrt(3)*sin(zeta))*sin(alpha - phi)));

    By_NV = @(Bz,phi) 1/6*(sqrt(2)*(bz+Bz)*sin(beta)*(2 - sqrt(3)*sin(zeta)) - ...
        Bperp*cos(alpha - phi)*(sqrt(3)*sin(zeta) + sin(beta)*(sqrt(2) + ...
        sqrt(6)*sin(zeta))) + cos(beta)*(2*(bz+Bz)*(1 - cos(zeta) + ...
        sqrt(3)*sin(zeta)) + Bperp*(2 + cos(zeta) - sqrt(3)*sin(zeta))* ...
        (cos(alpha - phi) - sin(alpha - phi))) + Bperp*(-sqrt(3)*sin(zeta) ...
        + sin(beta)*(sqrt(2) + sqrt(6)*sin(zeta)))*sin(alpha - phi) + ...
        cos(zeta)*(sqrt(2)*(bz+Bz)*sin(beta) + Bperp*cos(alpha - phi)*(-3 + ...
        sqrt(2)*sin(beta)) - Bperp*(3 + sqrt(2)*sin(beta))*sin(alpha - phi)));

    Bz_NV = @(Bz,phi) 1/6*(-2*sqrt(2)*(bz+Bz)*(-1 + cos(zeta))*sin(beta) - ...
        sqrt(2)*Bperp*(1 + 2*cos(zeta))*sin(beta)*(cos(alpha - phi) - ...
        sin(alpha - phi)) + 2*sqrt(3)*Bperp*sin(zeta)*(cos(alpha - phi) ...
        + sin(alpha - phi)) + 2*cos(beta)*((bz+Bz) + 2*(bz+Bz)*cos(zeta) + ...
        Bperp*cos(alpha - phi) - Bperp*cos(zeta)*cos(alpha - phi) + ...
        Bperp*(-1 + cos(zeta))*sin(alpha - phi)));


    PL = ones(length(Bz_list), length(phi_list));
    dPLdBz = zeros(length(Bz_list), length(phi_list));
    dPLdphi = zeros(length(Bz_list), length(phi_list));

    % Evaluate PL
    for i = 1:length(Bz_list)
        Bz = Bz_list(i);
        for j = 1:length(phi_list)
            phi = phi_list(j);
            bx = Bx_NV(Bz,phi);
            by = By_NV(Bz,phi);
            bz_loc = Bz_NV(Bz,phi); % local variable to avoid shadowing input `bz`
            

            % Finite differences for field-component derivatives
            dbxdBz = (Bx_NV(Bz+dBz,phi)-bx) / dBz;
            dbydBz = (By_NV(Bz+dBz,phi)-by) / dBz; 
            dbzdBz = (Bz_NV(Bz+dBz,phi)-bz_loc) / dBz;

            dbxdphi = (Bx_NV(Bz,phi+dphi)-bx) / dphi;
            dbydphi = (By_NV(Bz,phi+dphi)-by) / dphi; 
            dbzdphi = (Bz_NV(Bz,phi+dphi)-bz_loc) / dphi;

            % PL model
            PL(i,j) = PL(i,j) - contrast / (1+((bx-bz_loc)/width)^2);
            PL(i,j) = PL(i,j) - contrast / (1+((bx-by)/width)^2);
            PL(i,j) = PL(i,j) - contrast / (1+((by-bz_loc)/width)^2);
            PL(i,j) = PL(i,j) - contrast / (1+((bx+bz_loc)/width)^2);
            PL(i,j) = PL(i,j) - contrast / (1+((bx+by)/width)^2);
            PL(i,j) = PL(i,j) - contrast / (1+((by+bz_loc)/width)^2);
            PL(i,j) = PL(i,j) - 2*contrast / (1+(bx/width)^2);
            PL(i,j) = PL(i,j) - 2*contrast / (1+(by/width)^2);
            PL(i,j) = PL(i,j) - 2*contrast / (1+(bz_loc/width)^2);

        
            % Derivatives dPL / dBz
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((bx-bz_loc)/width)^2).^2.*2*(bx-bz_loc)/width^2*(dbxdBz-dbzdBz);
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((bx-by)/width)^2).^2.*2*(bx-by)/width^2*(dbxdBz-dbydBz);
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((by-bz_loc)/width)^2).^2.*2*(by-bz_loc)/width^2*(dbydBz-dbzdBz);
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((bx+bz_loc)/width)^2).^2.*2*(bx+bz_loc)/width^2*(dbxdBz+dbzdBz);
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((bx+by)/width)^2).^2.*2*(bx+by)/width^2*(dbxdBz+dbydBz);
            dPLdBz(i,j) = dPLdBz(i,j) + contrast / (1+((by+bz_loc)/width)^2).^2.*2*(by+bz_loc)/width^2*(dbydBz+dbzdBz);
            dPLdBz(i,j) = dPLdBz(i,j) + 2*contrast / (1+(bx/width)^2).^2.*2*(bx)/width^2*(dbxdBz);
            dPLdBz(i,j) = dPLdBz(i,j) + 2*contrast / (1+(by/width)^2).^2.*2*(by)/width^2*(dbydBz);
            dPLdBz(i,j) = dPLdBz(i,j) + 2*contrast / (1+(bz_loc/width)^2).^2.*2*(bz_loc)/width^2*(dbzdBz);

            
            % Derivatives dPL/ dphi
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((bx-bz_loc)/width)^2).^2.*2*(bx-bz_loc)/width^2*(dbxdphi-dbzdphi);
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((bx-by)/width)^2).^2.*2*(bx-by)/width^2*(dbxdphi-dbydphi);
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((by-bz_loc)/width)^2).^2.*2*(by-bz_loc)/width^2*(dbydphi-dbzdphi);
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((bx+bz_loc)/width)^2).^2.*2*(bx+bz_loc)/width^2*(dbxdphi+dbzdphi);
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((bx+by)/width)^2).^2.*2*(bx+by)/width^2*(dbxdphi+dbydphi);
            dPLdphi(i,j) = dPLdphi(i,j) + contrast / (1+((by+bz_loc)/width)^2).^2.*2*(by+bz_loc)/width^2*(dbydphi+dbzdphi);
            dPLdphi(i,j) = dPLdphi(i,j) + 2*contrast / (1+(bx/width)^2).^2.*2*(bx)/width^2*(dbxdphi);
            dPLdphi(i,j) = dPLdphi(i,j) + 2*contrast / (1+(by/width)^2).^2.*2*(by)/width^2*(dbydphi);
            dPLdphi(i,j) = dPLdphi(i,j) + 2*contrast / (1+(bz_loc/width)^2).^2.*2*(bz_loc)/width^2*(dbzdphi);
        end
    end


    % Total variance including propagated sweep uncertainties
    sigma2 = sigma^2 + dPLdBz.^2 * sigma_Bz^2 + dPLdphi.^2*sigma_phi^2;
    
    % Model vs experimental - Gaussian log-likelihood
    S_mod = PL;
    s = size(S_mod);
    if s(1) ~= s(2) && all(s == size(S_exp'))
        S_mod = S_mod';
        sigma2 = sigma2';
    end

    % Prior distributions are assumed to be flat (no information about the
    % parameters), so the likelihood directly relates to the posterior
    % probability
    logL = -0.5 * sum((S_exp(:) - S_mod(:)).^2 ./ sigma2(:));
end