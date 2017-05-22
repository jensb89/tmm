classdef TMM < handle
    %TMM Transfer Matrix Method code for Matlab
    %   This code is a port of the Python Repository code from Steven
    %   Byrnes. Some functions of the original python code are maybe
    %   missing. Check out his repository: https://github.com/sbyrnes321/tmm
    %   Usage example: tmm = TMM()
    %                  [T,R]=tmm.coh_tmm('p',n_list, d_list, 0, lambda)
    %
    %   For information see the docstring of each function, and also see
    %   manual.pdf (should be included with the distribution, otherwise get it
    %   at http://sjbyrnes.com/fresnel_manual.pdf ). Physics background,
    %   conventions, and derivations are at https://arxiv.org/abs/1603.02720
    %   The most two important functions are:
    %   coh_tmm(...) -- the transfer-matrix-method calculation in the coherent
    %   case (i.e. thin films)
    %   inc_tmm(...) -- the transfer-matrix-method calculation in the incoherent
    %   case (i.e. films tens or hundreds of wavelengths thick, or whose
    %   thickness is not very uniform.) --- currently not implemented in
    %   Matlab
    
    
    properties
        EPSILON = eps %typical floating-point calculation error
    end
    
    methods
        function this=TMM()
        end
        
        function answer = is_forward_angle(this,n, theta)
            % if a wave is traveling at angle theta from normal in a medium with index n,
            % calculate whether or not this is the forward-traveling wave (i.e., the one
            % going from front to back of the stack, like the incoming or outgoing waves,
            % but unlike the reflected wave). For real n & theta, the criterion is simply
            % -pi/2 < theta < pi/2, but for complex n & theta, it's more complicated.
            % See https://arxiv.org/abs/1603.02720 appendix D. If theta is the forward
            % angle, then (pi-theta) is the backward angle and vice-versa.

            assert(real(n) * imag(n) >= 0,"For materials with gain, it's ambiguous which beam is incoming vs outgoing. See https://arxiv.org/abs/1603.02720 Appendix C.\n n: " + string(n) + "   angle: " + string(theta))
            ncostheta = n * cos(theta);
            if abs(imag(ncostheta)) > 100 * this.EPSILON
                % Either evanescent decay or lossy medium. Either way, the one that
                % decays is the forward-moving wave
                answer = (imag(ncostheta) > 0);
            else
                % Forward is the one with positive Poynting vector
                % Poynting vector is Re[n cos(theta)] for s-polarization or
                % Re[n cos(theta*)] for p-polarization, but it turns out they're consistent
                % so I'll just assume s then check both below
                answer = (real(ncostheta) > 0);
            end

            % double-check the answer ... can't be too careful!
            error_string = "It's not clear which beam is incoming vs outgoing. Weird index maybe?\n n: " + string(n) + "   angle: " + string(theta);
            if answer == true
                assert(imag(ncostheta) > -100 * this.EPSILON, error_string);
                assert(real(ncostheta) > -100 * this.EPSILON, error_string);
                assert(real(n * cos(conj(theta))) > -100 * this.EPSILON, error_string);
            else
                assert(ncostheta.imag < 100 * this.EPSILON, error_string)
                assert(ncostheta.real < 100 * this.EPSILON, error_string)
                assert(real(n * cos(theta.conjugate())) < 100 * this.EPSILON, error_string)
            end
        end
        
        function ret = snell(this,n_1, n_2, th_1)
            %return angle theta in layer 2 with refractive index n_2, assuming
            %it has angle th_1 in layer with refractive index n_1. Use Snell's law. Note
            %that "angles" may be complex!!
        
            % Important that the arcsin here is scipy.arcsin, not numpy.arcsin! (They
            % give different results e.g. for arcsin(2).)
            th_2_guess = asin(n_1*np.sin(th_1) / n_2);
            if this.is_forward_angle(n_2, th_2_guess)
                ret = th_2_guess;
            else
                ret =  pi - th_2_guess;
            end
        end
        
        function angles = list_snell(this, n_list, th_0)
            %return list of angle theta in each layer based on angle th_0 in layer 0,
            %using Snell's law. n_list is index of refraction of each layer. Note that
            %"angles" may be complex!!

            angles = asin(n_list(1)*sin(th_0) ./ n_list);
            % The first and last entry need to be the forward angle (the intermediate
            % layers don't matter, see https://arxiv.org/abs/1603.02720 Section 5)
            if ~this.is_forward_angle(n_list(1), angles(1))
                angles(1) = pi - angles(1);
            end
            if ~this.is_forward_angle(n_list(end), angles(end))
                angles(end) = pi - angles(end);
            end
        end
        
        function ret = interface_r(this,polarization, n_i, n_f, th_i, th_f)
            %%reflection amplitude (from Fresnel equations)
            %polarization is either "s" or "p" for polarization
            %n_i, n_f are (complex) refractive index for incident and final
            %th_i, th_f are (complex) propegation angle for incident and final
            %(in radians, where 0=normal). "th" stands for "theta".

            if polarization == 's'
                ret = ((n_i * cos(th_i) - n_f * cos(th_f)) / (n_i * cos(th_i) + n_f * cos(th_f)));
            elseif polarization == 'p'
                ret = ((n_f * cos(th_i) - n_i * cos(th_f)) / (n_f * cos(th_i) + n_i * cos(th_f)));
            else
                error("Polarization must be 's' or 'p'");
            end
        end
        
        
        function ret = interface_t(this,polarization, n_i, n_f, th_i, th_f)
            %transmission amplitude (frem Fresnel equations)
            %polarization is either "s" or "p" for polarization
            %n_i, n_f are (complex) refractive index for incident and final
            %th_i, th_f are (complex) propegation angle for incident and final
            %(in radians, where 0=normal). "th" stands for "theta".

            if polarization == 's'
                ret= 2 * n_i * cos(th_i) / (n_i * cos(th_i) + n_f * cos(th_f));
            elseif polarization == 'p'
                ret= 2 * n_i * cos(th_i) / (n_f * cos(th_i) + n_i * cos(th_f));
            else
                error("Polarization must be 's' or 'p'");
            end
        end

        function ret= R_from_r(this,r)
            %Calculate reflected power R, starting with reflection amplitude r.
            
            ret = abs(r).^2;
        end

        function ret = T_from_t(this,pol, t, n_i, n_f, th_i, th_f)
            %Calculate transmitted power T, starting with transmission amplitude t.
            %n_i,n_f are refractive indices of incident and final medium.
            %th_i, th_f are (complex) propegation angles through incident & final medium
            %(in radians, where 0=normal). "th" stands for "theta".
            %In the case that n_i,n_f,th_i,th_f are real, formulas simplify to
            %T=|t|^2 * (n_f cos(th_f)) / (n_i cos(th_i)).
            %See manual for discussion of formulas

            if pol == 's'
                ret = abs(t^2) * ((real(n_f*cos(th_f))) / real(n_i*cos(th_i)));
            elseif pol == 'p'
                ret= abs(t^2) * ((real(n_f*conj(cos(th_f)))) / real(n_i*conj(cos(th_i))));
            else
                error("Polarization must be 's' or 'p'");
            end
        end

        function ret = power_entering_from_r(this,pol, r, n_i, th_i)
            %Calculate the power entering the first interface of the stack, starting with
            %reflection amplitude r. Normally this equals 1-R, but in the unusual case
            %that n_i is not real, it can be a bit different than 1-R. See manual.
            %n_i is refractive index of incident medium.
            %th_i is (complex) propegation angle through incident medium
            %(in radians, where 0=normal). "th" stands for "theta".

            if pol == 's'
                ret = (real(n_i*cos(th_i)*(1+conj(r))*(1-r)) / real(n_i*cos(th_i)));
            elseif pol == 'p'
                ret = (real(n_i*conj(cos(th_i))*(1+r)*(1-conj(r))) / real(n_i*conj(cos(th_i))));
            else
                error("Polarization must be 's' or 'p'");
            end
        end
        
        function ret = interface_R(this,polarization, n_i, n_f, th_i, th_f)
            %Fraction of light intensity reflected at an interface.

            r = this.interface_r(polarization, n_i, n_f, th_i, th_f);
            ret =  this.R_from_r(r);
        end

        function ret = interface_T(this,polarization, n_i, n_f, th_i, th_f)
            %Fraction of light intensity transmitted at an interface.

            t = this.interface_t(polarization, n_i, n_f, th_i, th_f);
            ret = this.T_from_t(polarization, t, n_i, n_f, th_i, th_f);
        end
                
        function [T,R,t,r,power_entering,vw_list,kz_list,th_list] = coh_tmm(this,pol, n_list, d_list, th_0, lam_vac)
            % Main "coherent transfer matrix method" calc. Given parameters of a stack,
            % calculates everything you could ever want to know about how light
            % propagates in it. (If performance is an issue, you can delete some of the
            % calculations without affecting the rest.)
            % pol is light polarization, "s" or "p".
            % n_list is the list of refractive indices, in the order that the light would
            % pass through them. The 0'th element of the list should be the semi-infinite
            % medium from which the light enters, the last element should be the semi-
            % infinite medium to which the light exits (if any exits).
            % th_0 is the angle of incidence: 0 for normal, pi/2 for glancing.
            % Remember, for a dissipative incoming medium (n_list[0] is not real), th_0
            % should be complex so that n0 sin(th0) is real (intensity is constant as
            % a function of lateral position).
            % d_list is the list of layer thicknesses (front to back). Should correspond
            % one-to-one with elements of n_list. First and last elements should be "inf".
            % lam_vac is vacuum wavelength of the light.
            % Outputs the following
            % * T--transmitted wave power (as fraction of incident)
            % * R--reflected wave power (as fraction of incident)
            % * r--reflection amplitude
            % * t--transmission amplitude
            % * power_entering--Power entering the first layer, usually (but not always)
            %   equal to 1-R (see manual).
            % * vw_list-- n'th element is [v_n,w_n], the forward- and backward-traveling
            %   amplitudes, respectively, in the n'th medium just after interface with
            %   (n-1)st medium.
            % * kz_list--normal component of complex angular wavenumber for
            %   forward-traveling wave in each layer.
            % * th_list--(complex) propagation angle (in radians) in each layer


            % Input tests
            if (length(lam_vac)>1 || length(th_0)>1)
                error('This function is not vectorized; you need to run one calculation at a time (1 wavelength, 1 angle, etc.)')
            end
                                 
            %if (length(n_list.ndim ~= 1) || (d_list.ndim ~= 1) || (n_list.size ~= d_list.size)
            if length(n_list) ~= length(d_list)
                error("Problem with n_list or d_list!")
            end
            
            assert(d_list(1) == d_list(end) && d_list(1) == inf, 'd_list must start and end with inf!')
            assert(abs(imag(n_list(1)*sin(th_0))) < 100*this.EPSILON, 'Error in n0 or th0!')
            assert(this.is_forward_angle(n_list(1), th_0), 'Error in n0 or th0!')
            num_layers = length(n_list);

            % th_list is a list with, for each layer, the angle that the light travels
            % through the layer. Computed with Snell's law. Note that the "angles" may be
            % complex!
            th_list = this.list_snell(n_list, th_0);

            % kz is the z-component of (complex) angular wavevector for forward-moving
            % wave. Positive imaginary part means decaying.
            kz_list = 2 * pi * n_list .* cos(th_list) ./ lam_vac;

            % delta is the total phase accrued by traveling through a given layer.
            % Ignore warning about inf multiplication
            %olderr = sp.seterr(invalid='ignore')
            delta = kz_list .* d_list;
            %sp.seterr(**olderr)

            % For a very opaque layer, reset delta to avoid divide-by-0 and similar
            % errors. The criterion imag(delta) > 35 corresponds to single-pass
            % transmission < 1e-30 --- small enough that the exact value doesn't
            % matter.
            for i = 2:num_layers-1
                if imag(delta(i)) > 35
                    delta(i) = real(delta(i)) + 35j;
                    %if 'opacity_warning' not in globals():
                    %    global opacity_warning
                    %    opacity_warning = True
                        display(["Warning: Layers that are almost perfectly opaque ",...
                              "are modified to be slightly transmissive, ",...
                              "allowing 1 photon in 10^30 to pass through. It's ",...
                              "for numerical stability. This warning will not ",...
                              "be shown again."]);
                end
            end

            % t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
            % respectively, coming from i, going to j. Only need to calculate this when
            % j=i+1. (2D array is overkill but helps avoid confusion.)
            t_list = zeros(num_layers, num_layers);
            r_list = zeros(num_layers, num_layers);
            for i = 1:num_layers-1
                t_list(i,i+1) = this.interface_t(pol, n_list(i), n_list(i+1),...
                                            th_list(i), th_list(i+1));
                r_list(i,i+1) = this.interface_r(pol, n_list(i), n_list(i+1),...
                                            th_list(i), th_list(i+1));
            end
            
            % At the interface between the (n-1)st and nth material, let v_n be the
            % amplitude of the wave on the nth side heading forwards (away from the
            % boundary), and let w_n be the amplitude on the nth side heading backwards
            % (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
            % M_list[n]. M_0 and M_{num_layers-1} are not defined.
            % My M is a bit different than Sernelius's, but Mtilde is the same.
            M_list = zeros(num_layers, 2, 2);
            for i =2:num_layers-1
                M_list(i,:,:) = (1/t_list(i,i+1)) * [exp(-1j*delta(i)), 0; 0, exp(1j*delta(i))] * [1, r_list(i,i+1); r_list(i,i+1), 1];
            end
            
            Mtilde = [1, 0; 0, 1];
            
            for i = 2: num_layers-1
                Mtilde = Mtilde * squeeze(M_list(i,:,:));
            end
            Mtilde = [1, r_list(1,2); r_list(1,2), 1]/t_list(1,2) *  Mtilde;

            % Net complex transmission and reflection amplitudes
            r = Mtilde(2,1)/Mtilde(1,1);
            t = 1/Mtilde(1,1);

            % vw_list[n] = [v_n, w_n]. v_0 and w_0 are undefined because the 0th medium
            % has no left interface.
            vw_list = zeros(num_layers, 2);
            vw = [t; 0];
            vw_list(end,:) = vw';
            for i =num_layers-1:-1:1
                vw = squeeze(M_list(i,:,:)) * vw;
                vw_list(i,:) = vw';
            end

            % Net transmitted and reflected power, as a proportion of the incoming light
            % power.
            R = this.R_from_r(r);
            T = this.T_from_t(pol, t, n_list(1), n_list(end), th_0, th_list(end));
            power_entering = this.power_entering_from_r(pol, r, n_list(1), th_0);
        end
    end
    
end

