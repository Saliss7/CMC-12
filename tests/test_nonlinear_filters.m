% test_nonlinear_filters.m — checks for est_ekf/est_ukf. 
% Owner: Pessoa C
%   matlab -batch "addpath('src','tests'); test_nonlinear_filters"

function test_nonlinear_filters()

fprintf('=== test_nonlinear_filters ===\n');
rng(0);
p = params();

% Real plant when available, else the local EOM copy (stub returns all zeros).
x_probe = [0.0; 0.3; 0.4; 0.5];
if norm(plant_cartpole(0, x_probe, 1.0, p)) < 1e-12
    fprintf('plant_cartpole.m is still a stub, using local_dynamics.\n');
    ref = @(t,x,u,pp) local_dynamics(x, u, pp);
else
    fprintf('plant_cartpole.m implemented, validating against it.\n');
    ref = @plant_cartpole;
end

npass = 0; nfail = 0;

% 1. Jacobian vs central difference
thetas = [0, pi/4, -pi/4, pi/2, -pi/2, pi];
u_test = 1.3;
max_err = 0;
for th = thetas
    x0 = [0.2; 0.5; th; 0.7];
    Fa = jacobian_f(x0, u_test, p);
    Fn = num_jacobian(@(x) ref(0, x, u_test, p), x0);
    max_err = max(max_err, max(abs(Fa(:) - Fn(:))));
end
tol = 1e-6;
if max_err < tol
    fprintf('[PASS] Jacobian finite-diff: max err %.2e < %.0e\n', max_err, tol);
    npass = npass + 1;
else
    fprintf('[FAIL] Jacobian finite-diff: max err %.2e >= %.0e\n', max_err, tol);
    nfail = nfail + 1;
end

% 2. Smoke test: finite, symmetric-PSD P, observed states tracked, no divergence.
cases = { [0;0;pi-0.3;0], 'hanging-osc';
          [0;0;2.5;0.5],   'large-swing' };
obs_tol = 0.05;
div_tol = 100;

for fi = 1:2
    switch fi
        case 1, est = @est_ekf; name = 'EKF';
        case 2, est = @est_ukf; name = 'UKF';
    end
    for ci = 1:size(cases,1)
        x_true = cases{ci,1};
        [ok, msg] = run_filter_smoke(est, ref, x_true, p, obs_tol, div_tol);
        if ok
            fprintf('[PASS] %s smoke (%s): %s\n', name, cases{ci,2}, msg);
            npass = npass + 1;
        else
            fprintf('[FAIL] %s smoke (%s): %s\n', name, cases{ci,2}, msg);
            nfail = nfail + 1;
        end
    end
end

% 3. Interface: dbg when asked, empty on a dropped sample, 2-output still ok
z  = [0.01; deg2rad(5)];
for est = {@est_ekf, @est_ukf}
    fn = est{1};
    [xh, P, dbg] = fn(p.x0hat, p.P0, z, 0, p);
    ok3 = isequal(size(xh),[4 1]) && isequal(size(P),[4 4]) && ...
          isstruct(dbg) && isequal(size(dbg.nu),[2 1]) && ...
          isequal(size(dbg.S),[2 2]);
    [~, ~, dbg2] = fn(p.x0hat, p.P0, [NaN; NaN], 0, p);
    ok3 = ok3 && isempty(dbg2);
    [xh2, P2] = fn(p.x0hat, p.P0, z, 0, p); %#ok<ASGLU>
    ok3 = ok3 && isequal(size(xh2),[4 1]);
    if ok3
        fprintf('[PASS] %s interface (dbg + 2-output)\n', func2str(fn));
        npass = npass + 1;
    else
        fprintf('[FAIL] %s interface\n', func2str(fn));
        nfail = nfail + 1;
    end
end

fprintf('--------------------------------------\n');
fprintf('Total: %d passed, %d failed\n', npass, nfail);
if nfail > 0
    error('test_nonlinear_filters: %d test(s) failed.', nfail);
end
fprintf('ALL TESTS PASSED\n');

end


function [ok, msg] = run_filter_smoke(est, ref, x_true, p, obs_tol, div_tol)
    N        = 300;
    xhat     = x_true;           % smoke test starts at truth
    P        = p.P0;
    u        = 0;
    err_obs  = zeros(N,1);
    err_full = zeros(N,1);
    for k = 1:N
        if mod(k,10) == 0        % drop every 10th sample
            z = [NaN; NaN];
        else
            z = [x_true(1); x_true(3)] + [p.sigma_x; p.sigma_theta].*randn(2,1);
        end

        [xhat, P] = est(xhat, P, z, u, p);

        if any(~isfinite(xhat)) || any(~isfinite(P(:)))
            ok = false; msg = sprintf('non-finite at step %d', k); return
        end
        if norm(P - P','fro') > 1e-9
            ok = false; msg = sprintf('P not symmetric at step %d', k); return
        end
        [~, notpd] = chol(P + 1e-12*eye(4));
        if notpd
            ok = false; msg = sprintf('P not PSD at step %d', k); return
        end

        err_obs(k)  = norm(xhat([1 3]) - x_true([1 3]));
        err_full(k) = norm(xhat - x_true);
        x_true = rk4_step(@(t,x) ref(t,x,u,p), 0, x_true, p.dt);
    end
    obs  = mean(err_obs(end-49:end));
    full = mean(err_full(end-49:end));
    if full > div_tol
        ok = false; msg = sprintf('diverged (full err %.1f)', full); return
    end
    if obs > obs_tol
        ok = false; msg = sprintf('observed-state err %.3f > tol %.3f', obs, obs_tol); return
    end
    ok  = true;
    msg = sprintf('obs err %.4f, full err %.3f', obs, full);
end


% Synthetic cart-pole EOM, KEEP IN SYNC with plant_cartpole.m.
function xdot = local_dynamics(x, u, p)
    M = p.M; m = p.m; L = p.L; g = p.g; b = p.b;
    xd  = x(2); th = x(3); thd = x(4);
    a = M + m;
    d = m*L^2/3 + m*L^2;      % I + m L^2
    k = m*L;
    c = cos(th); s = sin(th);
    det  = a*d - k^2*c^2;
    rhs1 = u - b*xd - k*thd^2*s;
    rhs2 = g*k*s;
    xddot  = (d*rhs1 + k*c*rhs2) / det;
    thddot = (a*rhs2 + k*c*rhs1) / det;
    xdot = [xd; xddot; thd; thddot];
end


function J = num_jacobian(f, x)
    n = numel(x);
    f0 = f(x); %#ok<NASGU>
    m = numel(f(x));
    J = zeros(m, n);
    h = 1e-6;
    for j = 1:n
        dx = zeros(n,1); dx(j) = h;
        J(:,j) = (f(x+dx) - f(x-dx)) / (2*h);
    end
end


function x_next = rk4_step(f, t, x, dt)
    k1 = f(t,        x);
    k2 = f(t + dt/2, x + dt/2 * k1);
    k3 = f(t + dt/2, x + dt/2 * k2);
    k4 = f(t + dt,   x + dt   * k3);
    x_next = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
