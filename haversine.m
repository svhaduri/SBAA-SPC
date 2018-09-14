function d = haversine (lat1, lon1, lat2, lon2)
    dlat = degtorad(lat2-lat1);
    dlon = degtorad(lon2-lon1);
    lat1 = degtorad(lat1);
    lat2 = degtorad(lat2);
    a = (sin(dlat./2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon./2)).^2;
    c = 2 .* asin(sqrt(a));
%     R = 6372.8; % km
    R = 6372.8 * 1000; % meter
    d = R*c;
end