function [ Interpolant, dxdD, dxdT, dxdY ] = interpolateDifferentiateEos...
    ( D, T, Y, D1D, T1D, Y1D, V3D, OS )

  Interpolant = zeros( size( D, 1 ), 1 );
  dxdD = zeros( size( D, 1 ), 1 );
  dxdT = zeros( size( D, 1 ), 1 );
  dxdY = zeros( size( D, 1 ), 1 );
  
  for i = 1 : size( D, 1 )

    iD = find( D1D < D(i), 1, 'last' );
    iT = find( T1D < T(i), 1, 'last' );
    iY = find( Y1D < Y(i), 1, 'last' );
    
    if ( size(iD,1)*size(iT,1)*size(iY,1) ~= 1 )
        disp( 'Error in interpolateEos: Outside table boundary.' )
        if ( size( iD, 1 ) == 0 )
            disp( 'Density out of boundary.' )
            disp('Replaced by the boundary.')
            iD = 1;
        end
        if ( size( iT, 1 ) == 0  )
            disp('Temperature out of boundary.')
            disp('Replaced by the boundary.')
            iT = 1;
        end
        if ( size( iY, 1 ) == 0 )
            disp( 'Electron fraction out of boundary.' )
            disp('Replaced by the boundary.')
            iY = 1;
        end
%         break
    end

    dD  = log10( D(i) / D1D(iD) ) / log10( D1D(iD+1) / D1D(iD) );
    ddD = 1.0d0 / ( D(i) * log10( D1D(iD+1) / D1D(iD) ) );
    
    dT = log10( T(i) / T1D(iT) ) / log10( T1D(iT+1) / T1D(iT) );
    ddT = 1.0d0 / ( T(i) * log10( T1D(iT+1) / T1D(iT) ) );
    
    dY = ( Y(i) - Y1D(iY) ) / ( Y1D(iY+1) - Y1D(iY) );
    ddY = log(10.0)/ ( Y1D(iY+1) - Y1D(iY) );
    
    
    p000 = log10( V3D( iD  , iT  , iY   ) + OS );
    p100 = log10( V3D( iD+1, iT  , iY   ) + OS );
    p010 = log10( V3D( iD  , iT+1, iY   ) + OS );
    p110 = log10( V3D( iD+1, iT+1, iY   ) + OS );
    p001 = log10( V3D( iD  , iT  , iY+1 ) + OS );
    p101 = log10( V3D( iD+1, iT  , iY+1 ) + OS );
    p011 = log10( V3D( iD  , iT+1, iY+1 ) + OS );
    p111 = log10( V3D( iD+1, iT+1, iY+1 ) + OS );
    
    Interpolant(i)...
      = 10.0^(...
          (1.0d0 - dY)...
          * (   (1.0d0 - dD) * (1.0d0 - dT) * p000...
              +          dD  * (1.0d0 - dT) * p100...
              + (1.0d0 - dD) *          dT  * p010...
              +          dD  *          dT  * p110 )...
          +         dY...
          * (   (1.0d0 - dD) * (1.0d0 - dT) * p001...
              +          dD  * (1.0d0 - dT) * p101...
              + (1.0d0 - dD) *          dT  * p011...
              +          dD  *          dT  * p111 ) )...
          - OS;
      
      dxdD(i) ...
          = ( (Interpolant(i) ) * ddD ...
            * ( (1.0 - dY) * ( (dT - 1.0) * p000   ...
                                    +  ( 1.0 - dT) * p100   ...
                                    -             dT  * p010   ...
                                    +             dT  * p110 ) ...
                         + dY * ( (dT - 1.0) * p001   ...
                                    +  ( 1.0 - dT) * p101   ...
                                    -             dT  * p011   ...
                                    +             dT  * p111 ) ) );

      
      dxdT(i)...
          = ( ( Interpolant(i) ) * ddT ...
              * ( (1.0 - dY ) * ( (dD - 1.0) * p000   ...
                                       -             dD  * p100   ...
                                       +  ( 1.0 - dD) * p010   ...
                                       +             dD  * p110 ) ...
                            + dY * ( (dD - 1.0) * p001   ...
                                       -             dD  * p101   ...
                                       +   (1.0 - dD) * p011   ...
                                       +             dD  * p111 ) ) );

      
      dxdY(i)...
          = ( ( Interpolant(i) ) * ddY ...
            * ( ( (dD - 1.0)) * (1.0 - dT) * p000   ...
                -            dD  * (1.0 - dT) * p100   ...
                - ( 1.0 - dD) *           dT  * p010   ...
                -            dD  *           dT  * p110   ...
                +  (1.0 - dD) * (1.0 - dT) * p001   ...
                +            dD  * (1.0 - dT) * p101   ...
                +  (1.0 - dD) *           dT  * p011   ...
                +            dD  *           dT  * p111 ) );

    
  end

end

