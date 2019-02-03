function [p,T] = hrtest(x,N)
% Hermans-Rasson test for non-uniformity of circular data.
%
%   [P,T] = HRTEST(X[,N]) computes Hermans-Rasson's test for the
%   sample of angles in X. If X is a matrix, the test is computed
%   independently on each column of X. The p-value is computed by simulation.
%   The optional argument N determines the number of simulated distributions
%   used to evaluate the p-value.
%
%   H0: the population is uniformly distributed around the circle.
%   H1: the population is *not* uniformly distributed around the circle.
%
%   No assumptions are made about the number of modes or symmetry of the
%   distribution.
%
%   Input:
%     X - sample of angles (in radians),
%     N - the number of simulations used to compute the p-value (default: 1000).
%
%   Output:
%     p - p-value of Hermans-Rasson's test,
%     T - the value of Hermans-Rasson's test statistic.
%
%   Notes:
%     1. Hermans-Rasson actually describe a family of tests. The test
%        implemented here corresponds to their test statistic with beta = 2.
%        This reportedly provides a good compromise when the data are unimodal.
%     2. The Hermans-Rasson test performs poorly with symmetric distributions.
%        In symmetric cases you're better off transforming the data to be
%        unimodal and employing Raleigh's test (see Lander et al. 2018).
%     3. The data are assumed to be continuously distributed.
%
%   Refs:
%     1. Hermans M, Rasson J. (1985). A new Sobolev test for uniformity on the
%        circle. Biometrika. 72:698–702. doi: 10.1093/biomet/72.3.698.
%     2. Landler L, Ruxton GD, Malkemper EP. (2018). Circular data in biology:
%        advice for effectively implementing statistical procedures. Behavioral
%        ecology and sociobiology. 72(8) 128. doi: 10.1007/s00265-018-2538-y.

% 2019-02-03 - Shaun L. Cloherty <s.cloherty@ieee.org>

narginchk(1,2);

if nargin < 2
  N = 1e3;
end

assert(numel(N) == 1 && N > 0,'N must be a positive scalar.')

% calculate the Hermans-Rasson test statistic
T = tstat(x);

% calculate p-value by simulation
n = size(x,1);

urnd = 2*pi*rand(n,N); % uniformly distributed on [0,2*pi]
uT = tstat(urnd);

p = sum(bsxfun(@lt,uT(:),T))./N;

end

function T = tstat(x)
  % Hermans-Rasson test statistic, beta = 2
  n = size(x,2);

  T = NaN([1,n]);
  for col = 1:n
    dx = bsxfun(@minus,x(:,col),x(:,col)');
    dx = abs(dx(:));

    T(col) = sum(pi - abs(pi - dx) + 2.895*abs(sin(dx)));
  end
end
