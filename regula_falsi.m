function [x, xs, k] = regula_falsi(f, a, b, tol, maxit)
    fa = f(a);
    fb = f(b);

    % prvi približek (sekanta)
    x_prev = a - fa*(b - a)/(fb - fa);
    xs = x_prev;

    for k = 2:maxit
        fx = f(x_prev);

        % izbira novega intervala
        if fa*fx < 0
            b = x_prev;
            fb = fx;
        else
            a = x_prev;
            fa = fx;
        end

        % nov približek
        x_new = a - fa*(b - a)/(fb - fa);
        xs(end+1) = x_new; %#ok<AGROW>

        % kriterij ustavitve
        if abs(x_new - x_prev) < tol
            x = x_new;
            return;
        end

        x_prev = x_new;
    end

    x = x_prev;
end
