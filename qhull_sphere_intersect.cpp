// g++ -O3 -o qhull_sphere_intersect qhull_sphere_intersect.cpp -lqhull_r -lqhullcpp
#include <iostream>
#include <vector>
#include <cmath>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/QhullPoint.h>
#include <fstream>
#include <map>

using namespace orgQhull;

struct parameters{
  double dx = 0.02;
  double dy = 0.02;
  double dz = 0.02;
  double radius = 0.4;
  double sphereX = 0.0;
  double sphereY = 0.0; 
  double sphereZ = 0.0;
};

void writeVTK(const std::string& filename, const std::vector<double>& points, const QhullFacetList& facets) ;
void writeVTK(const std::string& filename, const std::vector<double>& points, 
              const QhullFacetList& facets, const parameters& para) ;

// Function to check if a grid cell intersects with sphere
bool intersectsSphere(const double x, const double y, const double z, const parameters & p){
    // Check if any corner of the grid cell is within radius of sphere center
    double corners[8][3] = {
        {x, y, z},
        {x+p.dx, y, z},
        {x, y+p.dy, z},
        {x+p.dx, y+p.dy, z},
        {x, y, z+p.dz},
        {x+p.dx, y, z+p.dz},
        {x, y+p.dy, z+p.dz},
        {x+p.dx, y+p.dy, z+p.dz}
    };
    
    bool hasInside = false;
    bool hasOutside = false;
    
    for(int i = 0; i < 8; i++) {
        double dist = std::sqrt(
            std::pow(corners[i][0] - p.sphereX, 2) +
            std::pow(corners[i][1] - p.sphereY, 2) + 
            std::pow(corners[i][2] - p.sphereZ, 2)
        );
        if(dist <= p.radius) hasInside = true;
        else hasOutside = true;
        
        if(hasInside && hasOutside) return true;
    }
    
    return false;
}

int calc_area(const parameters para) {

    std::vector<std::vector<double>> intersectingCells;
    
    // Generate grid and check intersections
    int count = 0;
    for(double x = -1.0; x <= 1.0; x += para.dx) {
        for(double y = -1.0; y <= 1.0; y += para.dy) {
            for(double z = -1.0; z <= 1.0; z += para.dz) {
                if(intersectsSphere(x, y, z, para)) {
                    intersectingCells.push_back({x, y, z});
                    count++;
                }
            }
        }
    }
    //std::cout << "Number of intersecting cells: " << count << "\n";

    // Calculate analytical surface area
    double analyticalArea = 4.0 * M_PI * para.radius * para.radius;
    double discreteArea = 0.0;
    
    // Use Qhull to compute convex hull and its surface area
    if(!intersectingCells.empty()) {
        std::vector<double> points;
        for(const auto& cell : intersectingCells) {
            points.push_back(cell[0]);
            points.push_back(cell[1]);
            points.push_back(cell[2]);
        }
#ifdef DEBUG
        std::cout << "Number of Points: " << points.size()/3 << std::endl;
#endif
        
        Qhull qhull;
        qhull.runQhull("", 3, points.size()/3, points.data(), "Qt FA");  // Added FA option for area calculation
        
        // Sum up areas of all facets
        discreteArea = 0.0;
        const QhullFacetList facets = qhull.facetList();
#ifdef DEBUG
        std::cout << "Number of Facets: " << facets.size() << std::endl;
#endif

        // Write the VTK file
        writeVTK("sphere_intersection.vtk", points, facets, para);

        for(const QhullFacet& facet : facets) {
            // Only include facets that are not at infinity
            if (!facet.isGood()) continue;
            
            // Get the vertices of the facet
            QhullVertexSet vertices = facet.vertices();
            //std::cout << "number of vertices in this qhull: " << vertices.size() << std::endl;
            if(vertices.size() < 3) continue;
            
            // Calculate the area using spherical triangulation
            double facetArea = 0.0;
            auto it = vertices.begin();
            QhullVertex first = *it;
            QhullVertex second = *(++it);
            
            for(int i = 2; i < vertices.size(); i++) {
                QhullVertex third = *(++it);
                
                // Get the coordinates
                QhullPoint p1 = first.point();
                QhullPoint p2 = second.point();
                QhullPoint p3 = third.point();
                
                // Project points onto the unit sphere
                double len1 = std::sqrt(p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2]);
                double len2 = std::sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2]);
                double len3 = std::sqrt(p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2]);
                
                double u1[3] = {p1[0]/len1, p1[1]/len1, p1[2]/len1};
                double u2[3] = {p2[0]/len2, p2[1]/len2, p2[2]/len2};
                double u3[3] = {p3[0]/len3, p3[1]/len3, p3[2]/len3};
                
                // Calculate the angles of the spherical triangle
                double cos_a = u2[0]*u3[0] + u2[1]*u3[1] + u2[2]*u3[2];
                double cos_b = u1[0]*u3[0] + u1[1]*u3[1] + u1[2]*u3[2];
                double cos_c = u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
                
                // Clamp values to avoid numerical issues
                cos_a = std::max(-1.0, std::min(1.0, cos_a));
                cos_b = std::max(-1.0, std::min(1.0, cos_b));
                cos_c = std::max(-1.0, std::min(1.0, cos_c));
                
                double a = std::acos(cos_a);
                double b = std::acos(cos_b);
                double c = std::acos(cos_c);
                
                // Calculate the spherical excess (area of the spherical triangle on unit sphere)
                double s = (a + b + c) / 2.0;
                double spherical_excess = 4.0 * std::atan(std::sqrt(
                    std::tan(s/2) * std::tan((s-a)/2) * std::tan((s-b)/2) * std::tan((s-c)/2)
                ));
                
                // Scale by radius squared to get the actual area on our sphere
                facetArea += spherical_excess * para.radius * para.radius;
            }
            
            discreteArea += facetArea;
        }
        
#ifdef DEBUG
        std::cout << "Convex hull vertices: " << qhull.vertexCount() << std::endl;
#endif
    }

    // Calculate percentage error
    double error = std::abs(discreteArea - analyticalArea) / analyticalArea * 100.0;
    
#ifdef DEBUG
    std::cout << "Number of intersecting cells: " << count << std::endl;
    std::cout << "Discrete surface area: " << discreteArea << std::endl;
    std::cout << "Analytical surface area: " << analyticalArea << std::endl;
    std::cout << "Percentage error: " << error << "%" << std::endl;
#endif
    printf("%10.2f%10.2f%10.2f%10.2f%10.2f\n", para.radius , para.dx , analyticalArea , discreteArea , error);

    return 0;
}

/**
 * Writes the intersecting points and facets to a VTK file with grid cell indices
 * @param filename The name of the output VTK file
 * @param points Vector of point coordinates (x1,y1,z1,x2,y2,z2,...)
 * @param facets List of facets from Qhull
 * @param para Parameters including grid spacing
 */
void writeVTK(const std::string& filename, const std::vector<double>& points, 
              const QhullFacetList& facets, const parameters& para) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write VTK header
    outFile << "# vtk DataFile Version 3.0\n";
    outFile << "Qhull Sphere Intersection\n";
    outFile << "ASCII\n";
    outFile << "DATASET POLYDATA\n";

    // Write points
    int numPoints = points.size() / 3;
    outFile << "POINTS " << numPoints << " float\n";
    for (int i = 0; i < numPoints; i++) {
        outFile << points[i*3] << " " << points[i*3+1] << " " << points[i*3+2] << "\n";
    }

    // Count valid facets and their vertices for POLYGONS section
    int validFacets = 0;
    int totalVertices = 0;
    for (const QhullFacet& facet : facets) {
        if (!facet.isGood()) continue;
        
        QhullVertexSet vertices = facet.vertices();
        if (vertices.size() < 3) continue;
        
        validFacets++;
        totalVertices += vertices.size();
    }

    // Write polygons (facets)
    outFile << "POLYGONS " << validFacets << " " << (validFacets + totalVertices) << "\n";
    
    // Store facet indices for later use in cell data
    std::vector<std::vector<int>> facetGridIndices;
    
    // Write facets
    for (const QhullFacet& facet : facets) {
        if (!facet.isGood()) continue;
        
        QhullVertexSet vertices = facet.vertices();
        if (vertices.size() < 3) continue;
        
        outFile << vertices.size();
        std::vector<int> vertexIndices;
        
        for (const QhullVertex& vertex : vertices) {
            // Get the point coordinates from the vertex
            QhullPoint point = vertex.point();
            
            // Find the index of this vertex in our points array
            double coords[3] = {point[0], point[1], point[2]};
            
            // Find the matching index
            int vertexIndex = -1;
            for (int i = 0; i < numPoints; i++) {
                if (std::abs(points[i*3] - coords[0]) < 1e-6 &&
                    std::abs(points[i*3+1] - coords[1]) < 1e-6 &&
                    std::abs(points[i*3+2] - coords[2]) < 1e-6) {
                    vertexIndex = i;
                    break;
                }
            }
            
            if (vertexIndex == -1) {
                std::cerr << "Error: Could not find matching vertex index for point ("
                          << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;
                vertexIndex = 0; // Fallback
            }
            
            outFile << " " << vertexIndex;
            vertexIndices.push_back(vertexIndex);
        }
        outFile << "\n";
        
        // Store the grid indices for this facet
        facetGridIndices.push_back(vertexIndices);
    }
    
    // Write cell data (grid indices)
    outFile << "CELL_DATA " << validFacets << "\n";
    
    // Write i-indices
    outFile << "SCALARS i_index int 1\n";
    outFile << "LOOKUP_TABLE default\n";
    for (const auto& indices : facetGridIndices) {
        // Calculate average i-index for the facet
        double avgI = 0.0;
        for (int idx : indices) {
            double x = points[idx*3];
            int i = static_cast<int>((x + 1.0) / para.dx);
            avgI += i;
        }
        avgI /= indices.size();
        outFile << static_cast<int>(std::round(avgI)) << "\n";
    }
    
    // Write j-indices
    outFile << "SCALARS j_index int 1\n";
    outFile << "LOOKUP_TABLE default\n";
    for (const auto& indices : facetGridIndices) {
        // Calculate average j-index for the facet
        double avgJ = 0.0;
        for (int idx : indices) {
            double y = points[idx*3 + 1];
            int j = static_cast<int>((y + 1.0) / para.dy);
            avgJ += j;
        }
        avgJ /= indices.size();
        outFile << static_cast<int>(std::round(avgJ)) << "\n";
    }
    
    // Write k-indices
    outFile << "SCALARS k_index int 1\n";
    outFile << "LOOKUP_TABLE default\n";
    for (const auto& indices : facetGridIndices) {
        // Calculate average k-index for the facet
        double avgK = 0.0;
        for (int idx : indices) {
            double z = points[idx*3 + 2];
            int k = static_cast<int>((z + 1.0) / para.dz);
            avgK += k;
        }
        avgK /= indices.size();
        outFile << static_cast<int>(std::round(avgK)) << "\n";
    }
    
    outFile.close();
    std::cout << "VTK file written to: " << filename << std::endl;
}

int main(){
	parameters para;
    for(double r = 0.8; r <= 0.8; r += 0.05){
        para.dx = 0.02;
        para.dy = 0.02;
        para.dz = 0.02;
        para.radius = r;
        para.sphereX = 0.0;
        para.sphereY = 0.0;
        para.sphereZ = 0.0;
		calc_area(para);
	}
}
