function y = centerOfSums(x, membershipValueOfX)

    x = x(:);
    membershipValueOfX = membershipValueOfX(:);

    total_area = sum( membershipValueOfX );
    y = sum( membershipValueOfX .* x ) / total_area;

end