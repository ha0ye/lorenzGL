/******************************************************************/
/* LorenzGL_verHY.app                                             */
/* Lorenz Attractor and lagged embedding reconstruction simulator */
/* May 17 2012                                                    */
/* Hao Ye                                                         */
/******************************************************************/
This program generates data for the canonical Lorenz Attractor (sigma = 10, rho = 28, beta = 8/3) 
to demonstrate reconstruction using lagged coordinate-embeddings.

=========== Keyboard Controls ======================================
'esc', 'q', 'Q'		-	quit program

'1'					-	show original attractor in 3-dimensions
'2'					-	show time series (projections of the attractor onto coordinate axes)
'3'					-	show lags of a single time series
'4'					-	show reconstruction using lags of a single time series
'5'					-	show reconstruction along with original attractor
'6'					-	show univariate prediction
'7'					-	show univariate prediction with time series
'8'					-	show cross-mapping with original attractor
'9'					-	show cross-mapping with time series

'`'					-	toggle drawing of tracing lines and extra labels
't', 'T'			-	toggle display of Takens' Theorem
'm', 'M'			-	toggle manifold label
'l', 'L'			-   toggle drawing of full time series (viewing modes 2,3,7,9)
'k', 'K'			-	toggle coloring of manifolds by distance to current point

'x', 'X'			-	toggle projection to x-axis (viewing mode 1)
					-	use lags of x for embedding (viewing mode 3-9)
'y', 'Y'			-	toggle projection to y-axis (viewing mode 1)
					-	use lags of y for embedding (viewing mode 3-9)
'z', 'Z'			-	toggle projection to z-axis (viewing mode 1)
					-	use lags of z for embedding (viewing mode 3-9)
'n', 'N'			-	toggle variable to xmap to (viewing mode 8 & 9)

'/'					-	toggle split view of time series and attractor (viewing modes 2,3)
' '					-	pause animation
'r'					-	restart animation from first frame
'f'					-	toggle full screen
'['					-	slow down animation
']'					-	speed up animation

=========== Mouse Controls =========================================
left-click & drag	-	rotate attractor (viewing mode 1, 4-9)
middle-click & drag	-	zoom in & out (viewing mode 1, 4-9) [may not fully work]
right-click & drag	-	move attractor (viewing mode 1, 4-9) [may not fully work]