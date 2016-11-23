From 8f76e1ff53e85a2cb1fbe90afc2e956538ce5f55 Mon Sep 17 00:00:00 2001
From: Jürgen Knödlseder <jurgen.knodlseder@irap.omp.eu>
Date: Sat, 30 May 2015 01:19:15 +0200
Subject: [PATCH] Implement rejection method for Monte Carlo also for HealPix maps.

---
 src/model/GModelSpatialDiffuseMap.cpp | 109 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++---------------------------
 1 file changed, 82 insertions(+), 27 deletions(-)

diff --git a/src/model/GModelSpatialDiffuseMap.cpp b/src/model/GModelSpatialDiffuseMap.cpp
index 347fc9d..5342eb5 100644
--- a/src/model/GModelSpatialDiffuseMap.cpp
+++ b/src/model/GModelSpatialDiffuseMap.cpp
@@ -31,6 +31,7 @@
 #include "GException.hpp"
 #include "GTools.hpp"
 #include "GMath.hpp"
+#include "GHealpix.hpp"
 #include "GModelSpatialDiffuseMap.hpp"
 #include "GModelSpatialRegistry.hpp"
 
@@ -322,18 +323,18 @@ double GModelSpatialDiffuseMap::eval_gradients(const GPhoton& photon) const
  * @param[in,out] ran Random number generator.
  * @return Sky direction.
  *
- * Returns a random sky direction according to the intensity distribution of
- * the model sky map. It makes use of a cache array that contains the
- * normalized cumulative flux values of the skymap. Using a uniform random
- * number, this cache array is scanned using a bi-section method to determine
- * the skymap pixel for which the position should be returned. To avoid
- * binning problems, the exact position within the pixel is set by a uniform
- * random number generator. In case of a 2D map, this is done by randomizing
- * the skymap pixel values (neglecting thus pixel distortions). In case of a
- * HealPix map, this is done by selecing a random sky position within the
- * HealPix map.
- *
- * @todo: Improve HealPix randomization
+ * Draws a random sky direction from the intensity distribution of the model
+ * sky map.
+ *
+ * The method makes use of a cache array that contains the normalized
+ * cumulative flux values of the skymap. Using a uniform random number, this
+ * cache array is scanned using a bi-section method to determine the skymap
+ * pixel for which the position should be returned.
+ *
+ * Within that pixel, a rejection method is used to draw a sky direction
+ * that follows the intensity distribution that is obtained when calling the
+ * interpolation operator. This assures that even for coarse binning of the
+ * sky map the simulation corresponds to the model.
  ***************************************************************************/
 GSkyDir GModelSpatialDiffuseMap::mc(const GEnergy& energy,
                                     const GTime&   time,
@@ -384,27 +385,54 @@ GSkyDir GModelSpatialDiffuseMap::mc(const GEnergy& energy,
         } // endif: had a 2D pixel
 
         // ... otherwise convert pixel into sky direction and randomize
-        // position. We use here a kluge to compute the radius that contains
-        // a HealPix pixel by setting this radius to twice the half angle
-        // subtended by the solid angle of the pixel. We then use a while
-        // loop to randomize the position within that circle and to retain
-        // only directions that result in the same pixel index.
+        // position. We use a while loop to randomize the position within a
+        // circle that enclosed the pixel and retain only directions that
+        // result in the same pixel index and that are compatible with the
+        // density distribution.
         else {
+
+            // Get pointer on HealPix projection
+            const GHealpix* healpix = static_cast<const GHealpix*>(m_map.projection());
+
+            // Get enclosing radius
+            double radius = healpix->max_pixrad();
+
+            // Initialize pixel centre
             dir = m_map.pix2dir(pixel);
+
+            // Get randomized pixel
             GSkyDir randomized_dir;
-            int     randomized_index = -1;
-            double  radius           = 2.0 * std::acos(1.0 - m_map.solidangle(index)/gammalib::twopi);
-            double  cosrad           = std::cos(radius);
-            while (randomized_index != index) {
+            double  cosrad = std::cos(radius);
+            while (true) {
+
+                // Get randomized sky direction
                 randomized_dir = dir;
                 double theta   = std::acos(1.0 - ran.uniform() * (1.0 - cosrad)) * gammalib::rad2deg;
                 double phi     = 360.0 * ran.uniform();
                 randomized_dir.rotate_deg(phi, theta);
-                randomized_index = m_map.dir2inx(randomized_dir);
-//std::cout << radius << " " << randomized_index << " " << index << " " << randomized_dir << std::endl;
-            }
+
+                // Skip if we are not in the actual pixel
+                if (m_map.dir2inx(randomized_dir) != index) {
+                    continue;
+                }
+
+                // Get map value at that sky direction
+                double value = m_map(randomized_dir);
+
+                // Get uniform random number up to the maximum
+                double uniform = ran.uniform() * m_mc_max[index];
+
+                // Exit loop if we're not larger than the map value
+                if (uniform <= value) {
+                    break;
+                }
+                
+            } // endwhile
+
+            // Store randomize sky position
             dir = randomized_dir;
-        }
+            
+        } // endelse: we had a HealPix map
 
     } // endif: there were pixels in sky map
 
@@ -772,7 +800,6 @@ void GModelSpatialDiffuseMap::prepare_map(void)
         // zero intensity in the skymap. Invalid pixels are also filtered.
         double sum = 0.0;
         for (int i = 0; i < npix; ++i) {
-            //double flux = m_map(i) * m_map.solidangle(i);
             double flux = m_map.flux(i);
             if (flux < 0.0 ||
                 gammalib::is_notanumber(flux) ||
@@ -808,10 +835,38 @@ void GModelSpatialDiffuseMap::prepare_map(void)
         // If we have a HealPix map then set radius to 180 deg
         if (m_map.projection()->code() == "HPX") {
 
+            // Get pointer on HealPix projection
+            const GHealpix* healpix =
+                static_cast<const GHealpix*>(m_map.projection());
+
             // Set the map radius to full sky
             m_radius = 180.0;
 
-        }
+            // Compute maximum value that may occur from bilinear
+            // interpolation within this pixel and push this value on the
+            // stack. We do this by checking values of all neighbours.
+            for (int i = 0; i < npix; ++i) {
+
+                // Get neighbours
+                std::vector<int> neighbours = healpix->neighbours(i);
+
+                // Loop over neighbours
+                double max = m_map(i);
+                for (int k = 0; k < neighbours.size(); ++k) {
+                    if (neighbours[k] != -1) {
+                        double value = m_map(neighbours[k]);
+                        if (value > max) {
+                            max = value;
+                        }
+                    }
+                }
+
+                // Store maximum
+                m_mc_max.push_back(max);
+            
+            } // endfor: looped over pixels
+
+        } // endif: Healpix projection
 
         // ... otherwise compute map centre and radius
         else {
--
libgit2 0.21.4

