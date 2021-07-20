addpath("~/chebfun");

  g = chebfun(@(t) besselj(0,t),[0,100]);
  plot(g), ylim([-.5 1])