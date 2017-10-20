/*
* Software License Agreement (BSD License)
*
*  Point Cloud Library (PCL) - www.pointclouds.org
*  Copyright (c) 2009-2012, Willow Garage, Inc.
*  Copyright (c) 2014, Stichting Centrum Wiskunde en Informatica.
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of Willow Garage, Inc. nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*/
#ifndef POINT_CLOUD_CODECV2_IMPL_HPP
#define POINT_CLOUD_CODECV2_IMPL_HPP

// point cloud compression from PCL
#include <pcl/compression/entropy_range_coder.h>
#include <pcl/compression/impl/entropy_range_coder.hpp>
#include <pcl/cloud_codec_v2/point_cloud_codec_v2.h>
#include <pcl/compression/point_coding.h>
#include <pcl/compression/impl/octree_pointcloud_compression.hpp>

//includes to do the ICP procedures
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <Eigen/Geometry>
#include <pcl/cloud_codec_v2/impl/rigid_transform_coding_impl.hpp>

const bool PROPOSED = true;
const int BLOCK_SIZE = 256;
const int METHOD = 2; // Combination1 : 0, Combination2 : 1, Combination3 : 2
//const unsigned int LAMDA_INDEX = 0;
//const double THRESHOLD = 8; // 0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32
const double LAMDA[7] = { 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0 };
//vector<char> ReconColorDataVector;

struct ColorData {
	unsigned char B, G, R;
	unsigned int IndexCount;
	unsigned int DifferenceValueSet[BLOCK_SIZE];
	unsigned int TotalDifference;
	vector <unsigned int> IndexList;
};

/*
struct index_data {
	unsigned int index, count;
	double abs_difference;
	unsigned int data_r, data_g, data_b;
	double data_y, data_u, data_v;
};

struct color_data {
	bool pixel_indexing;
	unsigned int data_r, data_g, data_b;
	double data_y, data_u, data_v;
	unsigned int data_index;
};

bool compare_count(const struct index_data &a, const struct index_data &b) {
	if (a.count != b.count) {
		return a.count > b.count;
	} else {
		return a.abs_difference < b.abs_difference;
	}
}*/

namespace pcl{

	namespace io{

		/*!
		* \brief  encoding routing based on overriding pcl octree codec written by Julius Kammerl
		*        extra features include 
		*        - centroid coding
		*        - color coding based on jpeg
		*        - scalable bitstream (not yet implemented)
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param  cloud_arg
		* \param  compressed_tree_data_out_arg
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void OctreePointCloudCodecV2<
			PointT, LeafT, BranchT, OctreeT>::encodePointCloud (
			const PointCloudConstPtr &cloud_arg,
			std::ostream& compressed_tree_data_out_arg)
		{
			unsigned char recent_tree_depth =
				static_cast<unsigned char> (getTreeDepth ());

			// CWI addition to prevent crashes as in original cloud codec
			deleteCurrentBuffer();
			deleteTree();

			// initialize octree
			setInputCloud (cloud_arg);

			// CWI added, when encoding hte output variable stores the (simplified icloud)
			output_ = PointCloudPtr(new PointCloud());

			// add point to octree
			addPointsFromInputCloud ();

			// make sure cloud contains points
			if (leaf_count_>0) {

				// color field analysis
				cloud_with_color_ = false;
				std::vector<pcl::PCLPointField> fields;
				int rgba_index = -1;
				rgba_index = pcl::getFieldIndex (*input_, "rgb", fields);
				if (rgba_index == -1)
				{
					rgba_index = pcl::getFieldIndex (*input_, "rgba", fields);
				}
				if (rgba_index >= 0)
				{
					point_color_offset_ = static_cast<unsigned char> (fields[rgba_index].offset);
					cloud_with_color_ = true;
				}

				// apply encoding configuration
				cloud_with_color_ &= do_color_encoding_;

				// if octree depth changed, we enforce I-frame encoding
				i_frame_ |= (recent_tree_depth != getTreeDepth ());// | !(iFrameCounter%10);

				// enable I-frame rate
				if (i_frame_counter_++==i_frame_rate_)
				{
					i_frame_counter_ =0;
					i_frame_ = true;
				}

				// increase frameID
				frame_ID_++;

				// do octree encoding
				if (!do_voxel_grid_enDecoding_)
				{
					point_count_data_vector_.clear ();
					point_count_data_vector_.reserve (cloud_arg->points.size ());
				}

				// initialize color encoding (including new color coding based on jpeg)
				if(!color_coding_type_){
					color_coder_.initializeEncoding ();
					color_coder_.setPointCount (static_cast<unsigned int> (cloud_arg->points.size ()));
					color_coder_.setVoxelCount (static_cast<unsigned int> (leaf_count_));
				}else
				{ // new jpeg color coding
					jp_color_coder_.initializeEncoding ();
					jp_color_coder_.setPointCount (static_cast<unsigned int> (cloud_arg->points.size ()));
					jp_color_coder_.setVoxelCount (static_cast<unsigned int> (leaf_count_));
				}
				// initialize point encoding (including new centroid encoder)
				point_coder_.initializeEncoding ();
				centroid_coder_.initializeEncoding ();
				point_coder_.setPointCount (static_cast<unsigned int> (cloud_arg->points.size ()));
				centroid_coder_.initializeEncoding ();
				centroid_coder_.setPointCount( static_cast<unsigned int> (object_count_));

				// serialize octree
				if (i_frame_)
					// i-frame encoding - encode tree structure without referencing previous buffer
					serializeTree (binary_tree_data_vector_, false);
				else
					// p-frame encoding - XOR encoded tree structure
					serializeTree (binary_tree_data_vector_, true);


				// write frame header information to stream
				writeFrameHeader (compressed_tree_data_out_arg);

				// apply entropy coding to the content of all data vectors and send data to output stream
				entropyEncoding(compressed_tree_data_out_arg, compressed_tree_data_out_arg);

				// reset object count
				object_count_ = 0;

				if (b_show_statistics_) // todo update for codec v2
				{
					float bytes_per_XYZ = static_cast<float> (compressed_point_data_len_) / static_cast<float> (point_count_);
					float bytes_per_color = static_cast<float> (compressed_color_data_len_) / static_cast<float> (point_count_);

					PCL_INFO ("*** POINTCLOUD ENCODING ***\n");
					PCL_INFO ("Frame ID: %d\n", frame_ID_);
					if (i_frame_)
						PCL_INFO ("Encoding Frame: Intra frame\n");
					else
						PCL_INFO ("Encoding Frame: Prediction frame\n");
					PCL_INFO ("Number of encoded points: %ld\n", point_count_);
					PCL_INFO ("XYZ compression percentage: %f%%\n", bytes_per_XYZ / (3.0f * sizeof(float)) * 100.0f);
					PCL_INFO ("XYZ bytes per point: %f bytes\n", bytes_per_XYZ);
					PCL_INFO ("Color compression percentage: %f%%\n", bytes_per_color / (sizeof (int)) * 100.0f);
					PCL_INFO ("Color bytes per point: %f bytes\n", bytes_per_color);
					PCL_INFO ("Size of uncompressed point cloud: %f kBytes\n", static_cast<float> (point_count_) * (sizeof (int) + 3.0f  * sizeof (float)) / 1024);
					PCL_INFO ("Size of compressed point cloud: %d kBytes\n", (compressed_point_data_len_ + compressed_color_data_len_) / (1024));
					PCL_INFO ("Total bytes per point: %f\n", bytes_per_XYZ + bytes_per_color);
					PCL_INFO ("Total compression percentage: %f\n", (bytes_per_XYZ + bytes_per_color) / (sizeof (int) + 3.0f * sizeof(float)) * 100.0f);
					PCL_INFO ("Compression ratio: %f\n\n", static_cast<float> (sizeof (int) + 3.0f * sizeof (float)) / static_cast<float> (bytes_per_XYZ + bytes_per_color));
				}
			} else {
				if (b_show_statistics_)
					PCL_INFO ("Info: Dropping empty point cloud\n");
				deleteTree();
				i_frame_counter_ = 0;
				i_frame_ = true;
			}
		}

		/*!
		* \brief  decoding routing based on overriding pcl octree codec written by Julius Kammerl
		*        extra features include 
		*        - centroid coding
		*        - color coding based on jpeg
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param  compressed_tree_data_in_arg
		* \param  cloud_arg  decoded point cloud
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::decodePointCloud (
			std::istream& compressed_tree_data_in_arg,
			PointCloudPtr &cloud_arg)
		{

			// synchronize to frame header
			syncToHeader(compressed_tree_data_in_arg);

			// initialize octree
			switchBuffers ();

			// added to prevent crashes as happens with original cloud codec
			deleteCurrentBuffer();
			deleteTree();
			setOutputCloud (cloud_arg);

			// color field analysis
			cloud_with_color_ = false;
			std::vector<pcl::PCLPointField> fields;
			int rgba_index = -1;
			rgba_index = pcl::getFieldIndex (*output_, "rgb", fields);
			if (rgba_index == -1)
				rgba_index = pcl::getFieldIndex (*output_, "rgba", fields);
			if (rgba_index >= 0)
			{
				point_color_offset_ = static_cast<unsigned char> (fields[rgba_index].offset);
				cloud_with_color_ = true;
			}

			// read header from input stream
			readFrameHeader (compressed_tree_data_in_arg);

			// set the right grid pattern to the JPEG coder
			jp_color_coder_ = ColorCodingJPEG<PointT>(75,color_coding_type_);

			// decode data vectors from stream
			entropyDecoding(compressed_tree_data_in_arg, compressed_tree_data_in_arg);

			// initialize color and point encoding
			if(!color_coding_type_)
				color_coder_.initializeDecoding ();
			else
				jp_color_coder_.initializeDecoding ();

			point_coder_.initializeDecoding ();
			centroid_coder_.initializeDecoding();

			// initialize output cloud
			output_->points.clear ();
			output_->points.reserve (static_cast<std::size_t> (point_count_));

			if (i_frame_)
				// i-frame decoding - decode tree structure without referencing previous buffer
				deserializeTree (binary_tree_data_vector_, false);
			else
				// p-frame decoding - decode XOR encoded tree structure
				deserializeTree (binary_tree_data_vector_, true);

			// assign point cloud properties
			output_->height = 1;
			output_->width = static_cast<uint32_t> (cloud_arg->points.size ());
			output_->is_dense = false;

			//! todo update for cloud codecV2
			if (b_show_statistics_)
			{
				float bytes_per_XYZ = static_cast<float> (compressed_point_data_len_) / static_cast<float> (point_count_);
				float bytes_per_color = static_cast<float> (compressed_color_data_len_) / static_cast<float> (point_count_);

				PCL_INFO ("*** POINTCLOUDV2 DECODING ***\n");
				PCL_INFO ("Frame ID: %d\n", frame_ID_);
				if (i_frame_)
					PCL_INFO ("Encoding Frame: Intra frame\n");
				else
					PCL_INFO ("Encoding Frame: Prediction frame\n");
				PCL_INFO ("Number of encoded points: %ld\n", point_count_);
				PCL_INFO ("XYZ compression percentage: %f%%\n", bytes_per_XYZ / (3.0f * sizeof (float)) * 100.0f);
				PCL_INFO ("XYZ bytes per point: %f bytes\n", bytes_per_XYZ);
				PCL_INFO ("Color compression percentage: %f%%\n", bytes_per_color / (sizeof (int)) * 100.0f);
				PCL_INFO ("Color bytes per point: %f bytes\n", bytes_per_color);
				PCL_INFO ("Size of uncompressed point cloud: %f kBytes\n", static_cast<float> (point_count_) * (sizeof (int) + 3.0f * sizeof (float)) / 1024.0f);
				PCL_INFO ("Size of compressed point cloud: %f kBytes\n", static_cast<float> (compressed_point_data_len_ + compressed_color_data_len_) / 1024.0f);
				PCL_INFO ("Total bytes per point: %d bytes\n", static_cast<int> (bytes_per_XYZ + bytes_per_color));
				PCL_INFO ("Total compression percentage: %f%%\n", (bytes_per_XYZ + bytes_per_color) / (sizeof (int) + 3.0f * sizeof (float)) * 100.0f);
				PCL_INFO ("Compression ratio: %f\n\n", static_cast<float> (sizeof (int) + 3.0f * sizeof (float)) / static_cast<float> (bytes_per_XYZ + bytes_per_color));
			}
		}

		/*!
		* \brief  helper function to compute the delta frames
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param  pcloud_arg_in input point cloud to simplify
		* \param  out_cloud  output simplified point cloud
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT>::simplifyPCloud(const PointCloudConstPtr &pcloud_arg_in, 
			PointCloudPtr &out_cloud )
		{
			// this is the octree coding part of the predictive encoding
			OctreePointCloudCompression<PointT,LeafT,BranchT,OctreeT> *octree_simplifier = new OctreePointCloudCompression<PointT,LeafT,BranchT,OctreeT>
				(
				MANUAL_CONFIGURATION,
				false,
				point_resolution_,
				octree_resolution_,
				true,
				0,
				true,
				color_bit_resolution_ /*,0,do_voxel_centroid_enDecoding_*/
				);

			// use bounding box and obtain the octree grid
			octree_simplifier->defineBoundingBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
			octree_simplifier->setInputCloud (pcloud_arg_in);
			octree_simplifier->addPointsFromInputCloud ();

			/////////////////////////////////////////////////////////////////////////
			//! initialize output cloud
			//PointCloudPtr out_cloud(new PointCloud( octree_simplifier.leaf_count_, 1));
			out_cloud->width = (uint32_t) octree_simplifier->getLeafCount();
			out_cloud->height = 1;
			out_cloud->points.reserve(octree_simplifier->getLeafCount());  
			// cant get access to the number of leafs
			////////////////////////////////////////////////////////////////////////////

			///////////// compute the simplified cloud by iterating the octree ////////
			octree::OctreeLeafNodeIterator<OctreeT> it_ = octree_simplifier->leaf_begin();
			octree::OctreeLeafNodeIterator<OctreeT> it_end = octree_simplifier->leaf_end();

			for(int l_index =0;it_ !=it_end; it_++, l_index++)
			{
				// new point for the simplified cloud
				PointT l_new_point;

				//! centroid for storage
				std::vector<int>& point_indices = it_.getLeafContainer().getPointIndicesVector();

				// if centroid coding, store centroid, otherwise add 
				if(!do_voxel_centroid_enDecoding_)
				{
					octree_simplifier->genLeafNodeCenterFromOctreeKey(it_.getCurrentOctreeKey(),l_new_point);
				}
				else
				{
					Eigen::Vector4f cent;
					pcl::compute3DCentroid<PointT>(*pcloud_arg_in, point_indices, cent);

					l_new_point.x = cent[0];
					l_new_point.y = cent[1];
					l_new_point.z = cent[2];
				}
				long color_r=0;
				long color_g=0;
				long color_b=0;

				//! compute average color
				for(int i=0; i< point_indices.size();i++)
				{
					color_r+=pcloud_arg_in->points[point_indices[i]].r;
					color_g+=pcloud_arg_in->points[point_indices[i]].g;
					color_b+=pcloud_arg_in->points[point_indices[i]].b;
				}

				l_new_point.r = (char) (color_r / point_indices.size());
				l_new_point.g = (char) (color_g / point_indices.size());
				l_new_point.b = (char) (color_b / point_indices.size());

				out_cloud->points.push_back(l_new_point);
			}
			//////////////// done computing simplified cloud and octree structure //////////////////
			delete octree_simplifier;
		};
		/*!
		* \brief  helper function to  generate macroblock tree for computing shared macroblocks
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param  in_cloud input point cloud
		* \return pointer to a macroblock tree structure
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> OctreePointCloudCompression<PointT,LeafT,BranchT,OctreeT> *  
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT>::generate_macroblock_tree(PointCloudConstPtr  in_cloud)
		{
			MacroBlockTree *tree = new MacroBlockTree
				(
				MANUAL_CONFIGURATION,
				false,
				point_resolution_,
				octree_resolution_ * macroblock_size_,
				true,
				0,
				true,
				color_bit_resolution_
				/*,0,do_voxel_centroid_enDecoding_*/
				);

			// I frame coder
			tree->defineBoundingBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
			tree->setInputCloud (in_cloud); /* assume right argument was given, either the original p cloud or an i cloud */
			tree->addPointsFromInputCloud ();
			return tree;
		}

		/*!
		* \brief  helper function to  do ICP prediction between shared macroblocks
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param  i_cloud macroblock icloud
		* \param  p_cloud macroblock pcloud
		* \param  rigid_transform  output rigid transform
		* \param  has_converged  true if converged otherwise false
		* \param  rgb_offsets values for rgb offsets (if enabled)
		* pointer to a macroblock tree structure
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void 
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT>::do_icp_prediction(
			PointCloudPtr i_cloud,
			PointCloudPtr p_cloud,
			Eigen::Matrix4f &rigid_transform,
			bool & has_converged,
			char *rgb_offsets
			)
		{
			// check if it makes sense to do icp, i.e. input and output points approx. equal
			bool do_icp =  p_cloud->size() > 6 ? 
				(p_cloud->size() < i_cloud->size() * 2) 
				&&  (p_cloud->size() >= i_cloud->size() * 0.5) : false;

			if(!do_icp)
			{
				has_converged = false;
				return;
			}

			// color variance filter
			double in_av[3]={0,0,0};
			double out_av[3]={0,0,0};
			double in_var=0;
			double out_var=0;

			if(do_icp)
			{
				// create references to ease acccess via [] operator in the cloud
				pcl::PointCloud<PointT> & rcloud_out = *p_cloud;
				pcl::PointCloud<PointT> & rcloud_in = *i_cloud;

				for(int i=0; i<rcloud_in.size();i++)
				{
					in_av[0]+= (double) rcloud_in[i].r;
					in_av[1]+= (double) rcloud_in[i].g;
					in_av[2]+= (double) rcloud_in[i].b;
				}

				in_av[0]/=rcloud_in.size();
				in_av[1]/=rcloud_in.size();
				in_av[2]/=rcloud_in.size();

				// variance
				for(int i=0; i<rcloud_in.size();i++)
				{
					double val= (rcloud_in[i].r - in_av[0]) * (rcloud_in[i].r - in_av[0]) + 
						(rcloud_in[i].g - in_av[1]) * (rcloud_in[i].g - in_av[1]) +
						(rcloud_in[i].b - in_av[2]) * (rcloud_in[i].b - in_av[2]);

					in_var+=val;
				}
				in_var/=(3*rcloud_in.size());

				//
				for(int i=0; i<rcloud_out.size();i++)
				{
					out_av[0]+= (double) rcloud_out[i].r;
					out_av[1]+= (double) rcloud_out[i].g;
					out_av[2]+= (double) rcloud_out[i].b;
				}
				out_av[0]/=rcloud_out.size();
				out_av[1]/=rcloud_out.size();
				out_av[2]/=rcloud_out.size();

				for(int i=0; i<rcloud_out.size();i++)
				{
					double val= (rcloud_out[i].r - out_av[0]) * (rcloud_out[i].r - out_av[0]) + 
						(rcloud_out[i].g - out_av[1]) * (rcloud_out[i].g - out_av[1]) +
						(rcloud_out[i].b - out_av[2]) * (rcloud_out[i].b - out_av[2]);

					out_var+=val;
				}
				out_var/=(3*rcloud_out.size());

				// for segments with large variance, skip the icp prediction
				if(in_var > icp_var_threshold_ || out_var > icp_var_threshold_)
					do_icp = false;
			}


			char &rgb_offset_r=rgb_offsets[0]; 
			char &rgb_offset_g=rgb_offsets[1]; 
			char &rgb_offset_b=rgb_offsets[2]; 

			if(do_icp_color_offset_){
				if(std::abs(out_av[0] - in_av[0]) < 32)
					rgb_offset_r = (char)(out_av[0] - in_av[0]);
				if(std::abs(out_av[1] - in_av[1]) < 32)
					rgb_offset_g = (char)(out_av[1] - in_av[1]);
				if(std::abs(out_av[2] - in_av[2]) < 32)
					rgb_offset_b = (char)(out_av[2] - in_av[2]);
			}

			if(!do_icp)
			{
				has_converged = false;
				return;
			}

			// do the actual icp transformation
			pcl::IterativeClosestPoint<PointT, PointT> icp;

			icp.setInputSource(i_cloud);
			icp.setInputTarget(p_cloud);

			icp.setMaximumIterations (icp_max_iterations_);
			// Set the transformation epsilon (criterion 2)
			icp.setTransformationEpsilon (icp_max_iterations_);
			// Set the euclidean distance difference epsilon (criterion 3)
			icp.setEuclideanFitnessEpsilon (transformationepsilon_);

			pcl::PointCloud<PointT> Final;
			icp.align(Final);

			// compute the fitness for the colors

			if(icp.hasConverged() && icp.getFitnessScore() < point_resolution_ * 2)
			{
				has_converged = true;
				rigid_transform = icp.getFinalTransformation();
				return;
			}
			has_converged = false;
		}

		/** \brief generate a point cloud Delta to output stream
		* \param icloud_arg  point cloud to be used a I frame
		* \param pcloud_arg  point cloud to be encoded as a pframe
		* \param out_cloud_arg [out] the predicted frame 
		* \param i_coded_data intra encoded data 
		* \param p_coded_data inter encoded data
		* \param icp_on_original  (option to do icp on original or simplified clouds)
		* \param write_out_cloud  (flag to write the output cloud to out_cloud_arg)
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void 
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT>::generatePointCloudDeltaFrame( 
			const PointCloudConstPtr &icloud_arg /* icloud is already stored but can also be given as argument */,
			const PointCloudConstPtr &pcloud_arg,
			PointCloudPtr& out_cloud_arg, /* write the output cloud so we can assess the quality resulting from this algorithm */
			std::ostream& i_coded_data, 
			std::ostream& p_coded_data, 
			bool icp_on_original,
			bool write_out_cloud)
		{
			// intra coded points storage (points that cannot be predicted)
			typename pcl::PointCloud<PointT>::Ptr intra_coded_points(new pcl::PointCloud<PointT>());

			// initialize simplified predictive cloud
			PointCloudPtr simp_pcloud(new PointCloud());

			if(!icp_on_original){
				simplifyPCloud(pcloud_arg, simp_pcloud);
			}

			// keep track of the prediction statistics
			long macro_block_count=0;
			long shared_macroblock_count=0;
			long convergence_count=0;

			// initialize output cloud
			out_cloud_arg->height=1;
			out_cloud_arg->width =0;

			////////////// generate the octree for the macroblocks //////////////////////////////
			MacroBlockTree * i_block_tree = generate_macroblock_tree(icloud_arg);
			MacroBlockTree * p_block_tree = generate_macroblock_tree(icp_on_original ? pcloud_arg:simp_pcloud);

			//////////// iterate the predictive frame and find common macro blocks /////////////
			octree::OctreeLeafNodeIterator<OctreeT> it_predictive = p_block_tree->leaf_begin();
			octree::OctreeLeafNodeIterator<OctreeT> it_predictive_end = p_block_tree->leaf_end();

			for(;it_predictive!=it_predictive_end;++it_predictive)
			{
				macro_block_count++;
				const octree::OctreeKey current_key = it_predictive.getCurrentOctreeKey();
				pcl::octree::OctreeContainerPointIndices *i_leaf;
				typename pcl::PointCloud<PointT>::Ptr cloud_out (new pcl::PointCloud<PointT>(icp_on_original ? *pcloud_arg : *simp_pcloud , it_predictive.getLeafContainer().getPointIndicesVector()));

				if((i_leaf = i_block_tree->findLeaf(current_key.x,current_key.y,current_key.z)) != NULL)
				{
					typename pcl::PointCloud<PointT>::Ptr cloud_in (new pcl::PointCloud<PointT>(*icloud_arg, i_leaf->getPointIndicesVector()));
					shared_macroblock_count++;
					// shared block, do icp
					bool icp_success=false;
					Eigen::Matrix4f rt;
					char rgb_offsets[3]={0,0,0};

					do_icp_prediction(
						cloud_in,
						cloud_out,
						rt,
						icp_success,
						rgb_offsets
						);

					if(icp_success)
					{
						convergence_count++;


						// icp success, encode the rigid transform
						std::vector<int16_t> comp_dat;
						Eigen::Quaternion<float> l_quat_out;
						RigidTransformCoding<float>::compressRigidTransform(rt,comp_dat,l_quat_out);

						// write octree key, write rigid transform
						int16_t l_key_dat[3]={0,0,0};

						l_key_dat[0] = (int) current_key.x; 
						l_key_dat[1] = (int) current_key.y; 
						l_key_dat[2] = (int) current_key.z;

						// write the p coded data (we can add entropy encoding later)
						p_coded_data.write((const char *) l_key_dat ,3*sizeof(int16_t));
						p_coded_data.write((const char *) &comp_dat[0] ,comp_dat.size()*sizeof(int16_t));

						if(do_icp_color_offset_)
							p_coded_data.write((const char *) &rgb_offsets[0] ,3*sizeof(char));

						// following code is for generation of the predicted frame
						Eigen::Matrix4f mdec;
						Eigen::Quaternion<float> l_quat_out_dec;
						RigidTransformCoding<float>::deCompressRigidTransform(comp_dat, mdec,l_quat_out_dec);

						/* uncomment this code for debugging the rigid transform
						bool corrected_tf_matrix=false;
						for(int i=0; i<16;i++)
						{
						float diff=0;
						if((diff=std::abs(rt(i/4,i%4) - mdec(i/4,i%4))) > 0.01){
						//mdec(i/4,i%4) = rt(i/4,i%4);
						corrected_tf_matrix=true;
						std::cout << " error decoding rigid transform "
						<< comp_dat.size() << " index " << i/4 << "," 
						<< i%4 << std::endl;
						std::cout << " original " << rt << std::endl;
						std::cout << " decoded " <<  mdec << std::endl;
						std::cin.get();
						}
						}
						if(corrected_tf_matrix)
						std::cout << " matrix decoded not ok " <<std::endl;
						else
						std::cout << " matrix decoded ok " <<std::endl;
						*/
						if(write_out_cloud){
							// predicted point cloud
							pcl::PointCloud<PointT> manual_final;
							transformPointCloud<PointT, float>
								(  *cloud_in,
								manual_final, 
								mdec
								);

							// generate the output points
							for(int i=0; i < manual_final.size();i++){

								PointT &pt = manual_final[i];

								// color offset
								if(do_icp_color_offset_){
									pt.r+=rgb_offsets[0];
									pt.g+=rgb_offsets[1];
									pt.b+=rgb_offsets[2];
								}

								out_cloud_arg->push_back(pt);
							}
						}
					}
					else
					{
						// icp failed
						// add to intra coded points
						for(int i=0; i < cloud_out->size();i++){
							if(write_out_cloud)
								out_cloud_arg->push_back((*cloud_out)[i]);
							intra_coded_points->push_back((*cloud_out)[i]);
						}
					}
				}
				else
				{
					// exclusive block
					// add to intra coded points
					for(int i=0; i < cloud_out->size();i++)
					{
						if(write_out_cloud)
							out_cloud_arg->push_back((*cloud_out)[i]);
						intra_coded_points->push_back((*cloud_out)[i]);
					}
				}
			}  
			/* encode all the points that could not be predicted 
			from the previous coded frame in i frame manner with cluod_codec_v2 */
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT> intra_coder 
				(
				MANUAL_CONFIGURATION,
				false,
				point_resolution_,
				octree_resolution_,
				true,
				0,
				true,
				color_bit_resolution_,
				color_coding_type_, 
				do_voxel_centroid_enDecoding_
				);

			intra_coder.defineBoundingBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
			intra_coder.encodePointCloud(intra_coded_points, i_coded_data);

			//compute convergence statistics
			//compute the convergence statistics
			shared_macroblock_percentage_ =  (float)shared_macroblock_count / (float)macro_block_count ;
			shared_macroblock_convergence_percentage_ =  (float)convergence_count / (float)shared_macroblock_count;

			// clean up
			delete i_block_tree;
			delete p_block_tree;
		}

		/** \brief routine to encode a delta frame (predictive coding)
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param icloud_arg  point cloud to be used a I frame
		* \param pcloud_arg  point cloud to be encoded as a pframe
		* \param out_cloud_arg [out] the predicted frame 
		* \param i_coded_data intra encoded data 
		* \param p_coded_data inter encoded data
		* \param icp_on_original  (option to do icp on original or simplified clouds)
		* \param write_out_cloud  (flag to write the output cloud to out_cloud_arg)
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::encodePointCloudDeltaFrame (
			const PointCloudConstPtr &icloud_arg, const PointCloudConstPtr &pcloud_arg, PointCloudPtr &out_cloud_arg, 
			std::ostream& i_coded_data, std::ostream& p_coded_data, bool icp_on_original, bool write_out_cloud  )
		{
			// intra coded points storage (points that cannot be predicted)
			typename pcl::PointCloud<PointT>::Ptr intra_coded_points(new pcl::PointCloud<PointT>());

			// keep track of the prediction statistics
			long macro_block_count=0;
			long shared_macroblock_count=0;
			long convergence_count=0;

			// initialize predictive cloud
			PointCloudPtr simp_pcloud(new PointCloud());

			if(!icp_on_original){
				simplifyPCloud(pcloud_arg, simp_pcloud);
			}

			// initialize output cloud
			out_cloud_arg->height=1;
			out_cloud_arg->width =0;

			////////////// generate the octree for the macroblocks //////////////////////////////
			MacroBlockTree * i_block_tree = generate_macroblock_tree(icloud_arg);
			MacroBlockTree * p_block_tree = generate_macroblock_tree(icp_on_original ? pcloud_arg:simp_pcloud);

			//////////// iterate the predictive frame and find common macro blocks /////////////
			octree::OctreeLeafNodeIterator<OctreeT> it_predictive = p_block_tree->leaf_begin();
			octree::OctreeLeafNodeIterator<OctreeT> it_predictive_end = p_block_tree->leaf_end();

			for(;it_predictive!=it_predictive_end;++it_predictive)
			{
				const octree::OctreeKey current_key = it_predictive.getCurrentOctreeKey();
				pcl::octree::OctreeContainerPointIndices *i_leaf;
				typename pcl::PointCloud<PointT>::Ptr cloud_out (new pcl::PointCloud<PointT>(icp_on_original ? *pcloud_arg : *simp_pcloud , it_predictive.getLeafContainer().getPointIndicesVector()));
				macro_block_count++;

				if((i_leaf = i_block_tree->findLeaf(current_key.x,current_key.y,current_key.z)) != NULL)
				{
					shared_macroblock_count++;
					typename pcl::PointCloud<PointT>::Ptr cloud_in (new pcl::PointCloud<PointT>(*icloud_arg, i_leaf->getPointIndicesVector()));

					// shared block, do icp
					bool icp_success=false;
					Eigen::Matrix4f rt;
					char rgb_offsets[3]={0,0,0};

					do_icp_prediction(
						cloud_in,
						cloud_out,
						rt,
						icp_success,
						rgb_offsets
						);

					if(icp_success)
					{
						convergence_count++;
						// icp success, encode the rigid transform
						std::vector<int16_t> comp_dat;
						Eigen::Quaternion<float> l_quat_out;
						RigidTransformCoding<float>::compressRigidTransform(rt,comp_dat,l_quat_out);

						// write octree key, write rigid transform
						int16_t l_key_dat[3]={0,0,0};

						l_key_dat[0] = (int) current_key.x; 
						l_key_dat[1] = (int) current_key.y; 
						l_key_dat[2] = (int) current_key.z;

						// write the p coded data (we can add entropy encoding later)
						uint8_t chunk_size = (uint8_t) (3*sizeof(int16_t) + comp_dat.size()*sizeof(int16_t) + (do_icp_color_offset_ ? 3:0)); // size of the chunk
						p_coded_data.write((const char *) &chunk_size, sizeof(chunk_size));
						p_coded_data.write((const char *) l_key_dat, 3*sizeof(int16_t));
						p_coded_data.write((const char *) &comp_dat[0], comp_dat.size()*sizeof(int16_t));

						if(do_icp_color_offset_)
							p_coded_data.write((const char *) rgb_offsets ,3*sizeof(char));

						// following code is for generation of predicted frame
						Eigen::Matrix4f mdec;
						Eigen::Quaternion<float> l_quat_out_dec;
						RigidTransformCoding<float>::deCompressRigidTransform(comp_dat, mdec,l_quat_out_dec);

						// predicted point cloud
						pcl::PointCloud<PointT> manual_final;
						if(write_out_cloud){
							transformPointCloud<PointT, float>
								(  *cloud_in,
								manual_final, 
								mdec
								);
						}

						// generate the output points
						if(write_out_cloud){
							for(int i=0; i < manual_final.size();i++){

								PointT &pt = manual_final[i];

								// color offset
								if(do_icp_color_offset_){
									pt.r+=rgb_offsets[0];
									pt.g+=rgb_offsets[1];
									pt.b+=rgb_offsets[2];
								}

								out_cloud_arg->push_back(pt);
							}
						}
					}
					else
					{
						// icp failed
						// add to intra coded points
						for(int i=0; i < cloud_out->size();i++){
							if(write_out_cloud)
								out_cloud_arg->push_back((*cloud_out)[i]);
							intra_coded_points->push_back((*cloud_out)[i]);
						}
					}
				}
				else
				{
					// exclusive block
					// add to intra coded points
					for(int i=0; i < cloud_out->size();i++)
					{
						if(write_out_cloud)
							out_cloud_arg->push_back((*cloud_out)[i]);
						intra_coded_points->push_back((*cloud_out)[i]);
					}
				}
			}  

			/* encode all the points that could not be predicted 
			from the previous coded frame in i frame manner with cluod_codec_v2 */
			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT> intra_coder 
				(
				MANUAL_CONFIGURATION,
				false,
				point_resolution_,
				octree_resolution_,
				true,
				0,
				true,
				color_bit_resolution_,
				color_coding_type_, 
				do_voxel_centroid_enDecoding_
				);
			intra_coder.defineBoundingBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
			intra_coder.encodePointCloud(intra_coded_points, i_coded_data);

			//compute the convergence statistics
			shared_macroblock_percentage_ =  (float)shared_macroblock_count / (float)macro_block_count ;
			shared_macroblock_convergence_percentage_ =  (float)convergence_count / (float)shared_macroblock_count;

			// clean up
			delete i_block_tree;
			delete p_block_tree;
		}

		/** \brief routine to decode a delta frame (predictive coding)
		* \author Rufael Mekuria rufael.mekuria@cwi.nl
		* \param icloud_arg  point cloud to be used a I frame
		* \param cloud_out_arg  decoded frame
		* \param i_coded_data intra encoded data 
		* \param p_coded_data inter encoded data
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::decodePointCloudDeltaFrame( 
			const PointCloudConstPtr &icloud_arg, PointCloudPtr &cloud_out_arg, 
			std::istream& i_coded_data, std::istream& p_coded_data)
		{
			////////////// generate the macroblocks of the iframe //////////////////////////////
			MacroBlockTree * i_block_tree = generate_macroblock_tree(icloud_arg);

			///////////// inputs for decoding chunks of p data /////////////////////////////////
			int16_t keys_in[3]={0,0,0};
			uint8_t chunk_size=0;
			char rgb_offsets[3]={0,0,0};
			std::vector<int16_t> comp_dat_in;
			bool ended_unexpected=false;
			int decoded_mblocks=0;

			///////////////////////////decode p data ///////////////////////////
			while(p_coded_data.good())
			{
				// read the chunk size first
				p_coded_data.read((char *) &chunk_size, sizeof(uint8_t));

				if((chunk_size > 0) && p_coded_data.good()){

					// read the keys 
					p_coded_data.read((char *) keys_in, 3*sizeof(int16_t));

					// compute size of comp_dat (13 for vector mode a, 6 for quaternion mode)
					uint8_t comp_dat_size = chunk_size - 3*sizeof(int16_t) - (do_icp_color_offset_ ? 3:0);
					comp_dat_in.resize(comp_dat_size/sizeof(int16_t));

					//  load compdat
					p_coded_data.read((char *) &comp_dat_in[0] ,comp_dat_in.size()*sizeof(int16_t));

					//load color offsets
					if(do_icp_color_offset_)
						p_coded_data.read((char *) &rgb_offsets[0] ,3*sizeof(char));

					pcl::octree::OctreeContainerPointIndices *i_leaf;

					if((i_leaf = i_block_tree->findLeaf(keys_in[0],keys_in[1],keys_in[2])))
					{
						typename pcl::PointCloud<PointT>::Ptr cloud_in (new pcl::PointCloud<PointT>(*icloud_arg, i_leaf->getPointIndicesVector()));

						// following code is for generation of predicted frame
						Eigen::Matrix4f mdec;
						Eigen::Quaternion<float> l_quat_out_dec;
						RigidTransformCoding<float>::deCompressRigidTransform(comp_dat_in, mdec,l_quat_out_dec);

						// predicted point cloud
						pcl::PointCloud<PointT> p_chunk;

						// generate the pcoded data from the previous frame
						transformPointCloud<PointT, float>
							(  
							*cloud_in,
							p_chunk, 
							mdec
							);
						decoded_mblocks++;

						//decode the predictive points
						for(int j=0; j<p_chunk.size();j++)
						{
							PointT &p_point = (p_chunk)[j];

							if(do_icp_color_offset_)
							{
								p_point.r+=p_point.r + rgb_offsets[0];
								p_point.g+=p_point.g + rgb_offsets[1];
								p_point.b+=p_point.b + rgb_offsets[2];
							}

							cloud_out_arg->push_back(p_point);
						} 
					}
					else{
						std::cout << " failed decoding predictive frame, no corresponding i block " << std::endl;
					}
				}
				else
					break;
			}
			std::cout << " decoded: " << decoded_mblocks 
				<< " pblocks, resulting in " 
				<< cloud_out_arg->size() 
				<< " output points " 
				<< std::endl;

			// decode the intra coded points

			typename pcl::PointCloud<PointT>::Ptr intra_coded_points(new pcl::PointCloud<PointT>());

			OctreePointCloudCodecV2<PointT,LeafT,BranchT,OctreeT> intra_coder 
				(
				MANUAL_CONFIGURATION,
				false,
				point_resolution_,
				octree_resolution_,
				true,
				0,
				true,
				color_bit_resolution_,
				color_coding_type_, 
				do_voxel_centroid_enDecoding_
				);

			intra_coder.decodePointCloud(i_coded_data,intra_coded_points);

			// add the decoded input points to the output cloud
			for(int i=0; i< intra_coded_points->size();i++)
				cloud_out_arg->push_back((*intra_coded_points)[i]);

			// clean up
			delete i_block_tree;

		}

		// 
		struct local_color_enh{
			unsigned int r; 
			unsigned int g;
			unsigned int b; 
		};
		/*
		//! function for coding an enhancement layer, use the same input cloud!! 
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void 
		OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::encodeEnhancementLayer(const PointCloudConstPtr &cloud_arg, 
		std::ostream& compressed_tree_data_out_arg)
		{
		// get the leafs
		LeafNodeIterator leaf_it = leaf_begin();
		LeafNodeIterator leaf_end_it = leaf_end();

		StaticRangeCoder color_coder;
		StaticRangeCoder position_coder;

		// temporary storage of the occupancy codes
		std::vector<char> occupancy_codes;
		occupancy_codes.reserve(this->leaf_count_);

		// write the color code back
		std::vector<char> color_vector;

		// iterate the leafs
		while(leaf_it != leaf_end_it)
		{
		// get the point indices
		std::vector<int>& point_indices = leaf_it.getLeafContainer().getPointIndicesVector();

		// get the voxel grid size
		Eigen::Vector3f minpt;
		Eigen::Vector3f maxpt;
		getVoxelBounds(leaf_it, minpt,maxpt);

		// middle vertices
		Eigen::Vector3f midden = (maxpt - minpt) * 0.5 + minpt;

		// initialize the code position
		char code_pos=0u;

		// global and local color averages for the color
		unsigned int color_r=0;
		unsigned int color_g=0; 
		unsigned int color_b=0;

		// local color averages for the colors
		std::array<local_color_enh,8> colors_per_cell;
		std::array<int,8> occupancy_count={};

		// iterate all points and compute the occupancy code
		for(int i=0; i < point_indices.size(); i++)
		{
		unsigned char occ_code=0u;
		const PointT &l_point = (*cloud_arg)[point_indices[i]];

		// 
		if(l_point.x > midden.x())
		occ_code+=1;

		//
		if(l_point.y > midden.y())
		occ_code+=2;

		//
		if(l_point.z > midden.z())
		occ_code+=4;

		code_pos |= (1u << occ_code);

		// global color average //
		color_r+=l_point.r;
		color_g+=l_point.g;
		color_b+=l_point.b;

		// local color average //
		colors_per_cell[occ_code].r+= l_point.r;
		colors_per_cell[occ_code].g+= l_point.g;
		colors_per_cell[occ_code].b+= l_point.b;
		occupancy_count[occ_code]++;
		}

		// write the occupancy code back
		occupancy_codes.push_back( (char) code_pos);

		color_r/=point_indices.size();
		color_g/=point_indices.size();
		color_b/=point_indices.size();

		///////////  local color averages /////////////////////
		for(int k=0;k<8;k++)
		{
		if(occupancy_count[k])
		{
		colors_per_cell[k].r/= occupancy_count[k];
		colors_per_cell[k].g/= occupancy_count[k];
		colors_per_cell[k].b/= occupancy_count[k];

		int color_diff_r = (int) color_r - (int) colors_per_cell[k].r;

		if(color_diff_r < 128 && color_diff_r  >= -128)
		color_vector.push_back((char) color_diff_r);
		else
		color_vector.push_back(0);

		int color_diff_g = (int) color_g - (int) colors_per_cell[k].g;

		if(color_diff_g < 128 && color_diff_g  >= -128)
		color_vector.push_back((char) color_diff_g);
		else
		color_vector.push_back(0);

		int color_diff_b = (int) color_b - (int) colors_per_cell[k].b;

		if(color_diff_b < 128 && color_diff_b  >= -128)
		color_vector.push_back((char) color_diff_b);
		else
		color_vector.push_back(0);
		}
		}

		// increment 
		++leaf_it;
		}

		/// format of the enhancement layer  //////
		std::size_t pos_vec_size = occupancy_codes.size();
		std::size_t color_vec_size = color_vector.size();
		std::size_t comp_dat_size=0;

		compressed_tree_data_out_arg.write( (const char *) &pos_vec_size, sizeof(pos_vec_size ));
		compressed_tree_data_out_arg.write( (const char *) &color_vec_size , sizeof(color_vec_size));

		comp_dat_size+=position_coder.encodeCharVectorToStream(occupancy_codes, compressed_tree_data_out_arg);
		comp_dat_size+=color_coder.encodeCharVectorToStream(color_vector, compressed_tree_data_out_arg);
		}

		//! function for coding an enhancment layer
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
		OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::decodeEnhancementLayer(std::istream& compressed_tree_data_in_arg, PointCloudPtr &cloud_arg, PointCloudPtr &cloud_arg_enh)
		{      
		// decode all the occupancy codes
		std::size_t size_of_colors; 
		std::size_t size_of_occupancy_codes;

		std::vector<char> occ_codes; 
		std::vector<char> color_deltas; 

		compressed_tree_data_in_arg.read((char *) &size_of_colors, sizeof(size_of_colors));
		compressed_tree_data_in_arg.read((char *) &size_of_occupancy_codes,sizeof(size_of_occupancy_codes));

		occ_codes.resize(size_of_occupancy_codes);
		color_deltas.resize(size_of_colors);

		StaticRangeCoder pos_decoder;
		StaticRangeCoder color_decoder;

		pos_decoder.decodeStreamToCharVector(compressed_tree_data_in_arg, occ_codes);
		color_decoder.decodeStreamToCharVector(compressed_tree_data_in_arg, color_deltas);

		// get the leafs
		LeafNodeIterator leaf_it = leaf_begin();
		LeafNodeIterator leaf_end_it = leaf_end();

		std::vector<char>::iterator curr_occ = occ_codes.begin();
		std::vector<char>::iterator color_it= color_deltas.begin();

		assert(occ_codes.size() == this->leaf_count_);

		// apply all inverse operations
		// iterate the leafs
		while(leaf_it != leaf_end_it)
		{
		// do the decoding 
		OctreeKey ckey;
		ckey = leaf_it.getCurrentOctreeKey();

		for(unsigned char k=0;k<8;k++)
		{
		unsigned char ccode = (1u << k); 

		// test occupancy
		if( (*curr_code & ccode) == ccode){
		uint8_t x=0,y=0,z=0;
		PointT newPoint;

		int l = k;

		if(l >= 4)
		{
		//zpos
		z++;
		l-=4;
		}
		if(l>=2)
		{
		y++;
		l-=2; 
		}
		if(l>=1)
		{
		x++;
		}

		//
		newPoint.x = static_cast<float> ((static_cast<double> (ckey.x) + (x ? 0.75 : 0.25)) * resolution_ + min_x_);
		newPoint.y = static_cast<float> ((static_cast<double> (ckey.y) + (y ? 0.75 : 0.25)) * resolution_ + min_y_);
		newPoint.z = static_cast<float> ((static_cast<double> (ckey.z) + (z ? 0.75 : 0.25)) * resolution_ + min_z_);

		std::vector<int>& indicesV =leaf_it.getLeafContainer().getPointIndicesVector();
		int index_point = indicesV.size() > 0 ? indicesV[0] : -1;

		if(index_point != -1){
		newPoint.r = *color_it++ + (*cloud_arg[index_point]).r;
		newPoint.g = *color_it++ + (*cloud_arg[index_point]).g;
		newPoint.b = *color_it++ + (*cloud_arg[index_point]).b;
		}
		//enhancement layer code decoded
		cloud_arg_enh->push_back(newPoint);           
		}
		}

		// increment iterator
		++leaf_it;
		++curr_code;
		}
		}
		*/

		/////////////////// CODE OVERWRITING CODEC V1 to become codec V2 ////////////////////////////////

		//! write Frame header, extended for cloud codecV2
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::writeFrameHeader(std::ostream& compressed_tree_data_out_arg)
		{
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (frame_header_identifier_), strlen (frame_header_identifier_));
			//! use the base class and write some extended information on codecV2
			OctreePointCloudCompression<PointT, LeafT, BranchT, OctreeT>::writeFrameHeader(compressed_tree_data_out_arg);

			//! write additional fields for cloud codec v2
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&do_voxel_centroid_enDecoding_), sizeof (do_voxel_centroid_enDecoding_)); // centroid coding
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&do_connectivity_encoding_), sizeof (do_connectivity_encoding_));        // connectivity coding (not yet added)
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&create_scalable_bitstream_), sizeof (create_scalable_bitstream_));     // scalable bitstream
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&color_coding_type_), sizeof (color_coding_type_));                    // color coding type (added jpeg based modes)
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&macroblock_size_), sizeof (macroblock_size_));                       // settings for ICP based prediction
			compressed_tree_data_out_arg.write (reinterpret_cast<const char*> (&do_icp_color_offset_), sizeof (do_icp_color_offset_));              // settings for ICP based prediction
		};

		//! read Frame header, extended for cloud codecV2
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::readFrameHeader(std::istream& compressed_tree_data_in_arg)
		{
			//! use the base class and read some extended information on codecV2
			OctreePointCloudCompression<PointT>::readFrameHeader(compressed_tree_data_in_arg);

			//! read additional fields for cloud codec v2
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&do_voxel_centroid_enDecoding_), sizeof (do_voxel_centroid_enDecoding_));
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&do_connectivity_encoding_), sizeof (do_connectivity_encoding_));
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&create_scalable_bitstream_), sizeof (create_scalable_bitstream_));
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&color_coding_type_), sizeof (color_coding_type_));
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&macroblock_size_), sizeof (macroblock_size_));
			compressed_tree_data_in_arg.read (reinterpret_cast<char*> (&do_icp_color_offset_), sizeof (do_icp_color_offset_));
		};

		/** \brief Encode leaf node information during serialization, added jpeg color coding and voxel centroid coding
		* \param leaf_arg: reference to new leaf node
		* \param key_arg: octree key of new leaf node
		* \brief added centroids encoding compared to the original codec of PCL
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::serializeTreeCallback (LeafT &leaf_arg, const OctreeKey& key_arg)
		{
			// reference to point indices vector stored within octree leaf //
			const std::vector<int>& leafIdx = leaf_arg.getPointIndicesVector();
			double lowerVoxelCorner[3];
			PointT centroid;
			// //

			// calculate lower voxel corner based on octree key //
			lowerVoxelCorner[0] = static_cast<double> (key_arg.x) * resolution_ + min_x_;
			lowerVoxelCorner[1] = static_cast<double> (key_arg.y) * resolution_ + min_y_;
			lowerVoxelCorner[2] = static_cast<double> (key_arg.z) * resolution_ + min_z_;
			// //

			// points detail encoding
			if (!do_voxel_grid_enDecoding_)
			{
				// encode amount of points within voxel
				point_count_data_vector_.push_back (static_cast<int> (leafIdx.size ()));

				// differentially encode points to lower voxel corner
				point_coder_.encodePoints (leafIdx, lowerVoxelCorner, input_);

				if (cloud_with_color_) {
					// encode color of points
					if(!color_coding_type_) {
						color_coder_.encodePoints (leafIdx, point_color_offset_, input_);
					} else {
						jp_color_coder_.encodePoints (leafIdx, point_color_offset_, input_);
					}
				}
			}
			else // centroid or voxelgrid encoding
			{
				if (cloud_with_color_)
				{
					// encode average color of all points within voxel
					if(!color_coding_type_)
						color_coder_.encodeAverageOfPoints (leafIdx, point_color_offset_, input_);
					else
						jp_color_coder_.encodeAverageOfPoints (leafIdx, point_color_offset_, input_);
					std::vector<char> & l_colors = color_coding_type_ ? jp_color_coder_.getAverageDataVectorB() : color_coder_.getAverageDataVector();

					//! get the colors from the vector from the color coder, use to store to do temporal prediction
					centroid.r =  l_colors[l_colors.size() - 1];
					centroid.g =  l_colors[l_colors.size() - 2];
					centroid.b =  l_colors[l_colors.size() - 3];
				}
				if(!do_voxel_centroid_enDecoding_)
				{
					centroid.x = lowerVoxelCorner[0] + 0.5 * resolution_ ;
					centroid.y = lowerVoxelCorner[1] + 0.5 * resolution_ ;
					centroid.z = lowerVoxelCorner[2] + 0.5 * resolution_ ;
				}
				else
				{
					Eigen::Vector4f f_centroid;
					pcl::compute3DCentroid<PointT>(*(input_), leafIdx, f_centroid);

					centroid.x = f_centroid[0];
					centroid.y = f_centroid[1];
					centroid.z = f_centroid[2];

					centroid_coder_.encodePoint(lowerVoxelCorner, centroid);
				}
				// store the simplified cloud so that it can be used for predictive encoding
				output_->points.push_back(centroid);
			}
		}

		/** \brief Decode leaf nodes information during deserialization, added jpeg and voxel centroid coding
		* \param key_arg octree key of new leaf node
		* \brief added centroids encoding compared to the original codec
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::deserializeTreeCallback (LeafT&, const OctreeKey& key_arg)
		{
			double lowerVoxelCorner[3];
			std::size_t pointCount, i, cloudSize;
			PointT newPoint;
			pointCount = 1;

			if (!do_voxel_grid_enDecoding_)
			{
				// get current cloud size
				cloudSize = output_->points.size();

				// get amount of point to be decoded
				pointCount = *point_count_data_vector_iterator_;
				point_count_data_vector_iterator_++;

				// increase point cloud by amount of voxel points
				for (i = 0; i < pointCount; i++)
					output_->points.push_back (newPoint);

				// calculcate position of lower voxel corner
				lowerVoxelCorner[0] = static_cast<double> (key_arg.x) * resolution_ + min_x_;
				lowerVoxelCorner[1] = static_cast<double> (key_arg.y) * resolution_ + min_y_;
				lowerVoxelCorner[2] = static_cast<double> (key_arg.z) * resolution_ + min_z_;

				// decode differentially encoded points
				point_coder_.decodePoints (output_, lowerVoxelCorner, cloudSize, cloudSize + pointCount);
			}
			else
			{
				// decode the centroid or the voxel center
				if(do_voxel_centroid_enDecoding_)
				{
					PointT centroid_point;

					// calculcate position of lower voxel corner
					lowerVoxelCorner[0] = static_cast<double> (key_arg.x) * resolution_ + min_x_;
					lowerVoxelCorner[1] = static_cast<double> (key_arg.y) * resolution_ + min_y_;
					lowerVoxelCorner[2] = static_cast<double> (key_arg.z) * resolution_ + min_z_;

					centroid_coder_.decodePoint(newPoint,lowerVoxelCorner);
				}
				else
				{
					// calculate center of lower voxel corner
					newPoint.x = static_cast<float> ((static_cast<double> (key_arg.x) + 0.5) * resolution_ + min_x_);
					newPoint.y = static_cast<float> ((static_cast<double> (key_arg.y) + 0.5) * resolution_ + min_y_);
					newPoint.z = static_cast<float> ((static_cast<double> (key_arg.z) + 0.5) * resolution_ + min_z_);
				}

				// add point to point cloud
				output_->points.push_back (newPoint);
			}

			if (cloud_with_color_)
			{
				if (data_with_color_)
					// decode color information
					if(!color_coding_type_)
						color_coder_.decodePoints (output_, output_->points.size () - pointCount,
						output_->points.size (), point_color_offset_);
					else
						jp_color_coder_.decodePoints (output_, output_->points.size () - pointCount, output_->points.size (), point_color_offset_);
				else
					// set default color information

					color_coder_.setDefaultColor(output_, output_->points.size() - pointCount, output_->points.size(), point_color_offset_);
			}
		}

		/** \brief Synchronize to frame header
		* \param compressed_tree_data_in_arg: binary input stream
		* \brief use the new v2 header
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::syncToHeader (std::istream& compressed_tree_data_in_arg)
		{
			// sync to frame header
			unsigned int header_id_pos = 0;
			while (header_id_pos < strlen (frame_header_identifier_))
			{
				char readChar;
				compressed_tree_data_in_arg.read (static_cast<char*> (&readChar), sizeof (readChar));
				if (readChar != frame_header_identifier_[header_id_pos++])
					header_id_pos = (frame_header_identifier_[0]==readChar)?1:0;
			}
			//! read the original octree header
			OctreePointCloudCompression<PointT>::syncToHeader (compressed_tree_data_in_arg);
		};

		/** \brief Apply entropy encoding to encoded information and output to binary stream, added bitstream scalability and centroid encoding compared to  V1
		* \param compressed_tree_data_out_arg1 binary voxel output stream
		* \param compressed_tree_data_out_arg2 binary color output stream
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::entropyEncoding (std::ostream& compressed_tree_data_out_arg1, std::ostream& compressed_tree_data_out_arg2 )
		{
			uint64_t binary_tree_data_vector_size;
			uint64_t point_avg_color_data_vector_size;

			compressed_point_data_len_ = 0;
			compressed_color_data_len_ = 0;

			// encode binary octree structure
			binary_tree_data_vector_size = binary_tree_data_vector_.size ();
			compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&binary_tree_data_vector_size), sizeof (binary_tree_data_vector_size));
			compressed_point_data_len_ += entropy_coder_.encodeCharVectorToStream(binary_tree_data_vector_, compressed_tree_data_out_arg1);

			// store the number of bytes used for voxel encoding
			compression_performance_metrics[0] = compressed_point_data_len_ ;

			// encode centroids
			if(do_voxel_centroid_enDecoding_)
			{
				// encode differential centroid information
				std::vector<char>& point_diff_data_vector = centroid_coder_.getDifferentialDataVector ();
				uint32_t point_diff_data_vector_size = (uint32_t) point_diff_data_vector.size();
				compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&point_diff_data_vector_size), sizeof (uint32_t));
				compressed_point_data_len_ += entropy_coder_.encodeCharVectorToStream(point_diff_data_vector, compressed_tree_data_out_arg1);
			}

			// store the number of bytes used for voxel centroid encoding
			compression_performance_metrics[1] = compressed_point_data_len_ - compression_performance_metrics[0];

			// encode colors
			if (cloud_with_color_)
			{
				// encode averaged voxel color information
				std::vector<char>& pointAvgColorDataVector = color_coding_type_ ?  jp_color_coder_.getAverageDataVector () : color_coder_.getAverageDataVector ();

				if (PROPOSED) {
					if (METHOD == 0) {
						// encode Combination 1 begin
						vector<struct ColorData> ColorArray;
						vector<char> IndexDataStream, ColorDataStream;
						vector<char> RunVector, LengthVector;

						for (unsigned int idx = 0; idx < pointAvgColorDataVector.size(); idx += 3) {
							struct ColorData tempColorData;
							tempColorData.B = (unsigned char)pointAvgColorDataVector[idx];
							tempColorData.G = (unsigned char)pointAvgColorDataVector[idx + 1];
							tempColorData.R = (unsigned char)pointAvgColorDataVector[idx + 2];
							tempColorData.IndexCount = 1;
							tempColorData.TotalDifference = 0;
							ColorArray.push_back(tempColorData);

							if (idx + 3 == pointAvgColorDataVector.size() || ColorArray.size() == BLOCK_SIZE) {
								for (unsigned int i = 0; i < ColorArray.size() - 1; ++i) {
									for (unsigned int j = i + 1; j < ColorArray.size(); ++j) {
										int diff_b = (int)ColorArray[i].B - (int)ColorArray[j].B;
										int diff_g = (int)ColorArray[i].G - (int)ColorArray[j].G;
										int diff_r = (int)ColorArray[i].R - (int)ColorArray[j].R;

										double diff_rgb = pow((diff_r * diff_r) + (diff_g * diff_g) + (diff_b * diff_b), 0.5);
										if (diff_rgb <= THRESHOLD) {
											ColorArray[i].IndexList.push_back(j);
											ColorArray[j].IndexList.push_back(i);

											++ColorArray[i].IndexCount;
											++ColorArray[j].IndexCount;

											if (THRESHOLD != 0) {
												ColorArray[i].DifferenceValueSet[i] = 0;
												ColorArray[i].DifferenceValueSet[j] = (unsigned int)diff_rgb;
												ColorArray[j].DifferenceValueSet[j] = 0;
												ColorArray[j].DifferenceValueSet[i] = (unsigned int)diff_rgb;

												ColorArray[i].TotalDifference += (unsigned int)diff_rgb;
												ColorArray[j].TotalDifference += (unsigned int)diff_rgb;
											}
										}
									}
								}

								int BestIndex;
								unsigned char CurrentIndex = 0;
								unsigned int MaxSize, DifferenceValue;
								vector<char> tempIndexDataStream, tempColorDataStream;
								tempIndexDataStream.push_back(0);
								tempColorDataStream.resize(ColorArray.size());
								for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
									tempColorDataStream[i] = -1;
								}
								while (1) {
									BestIndex = -1;
									MaxSize = 0;
									DifferenceValue = 0;

									for (unsigned int i = 0; i < ColorArray.size(); ++i) {
										if (MaxSize < ColorArray[i].IndexCount || (THRESHOLD != 0 && MaxSize == ColorArray[i].IndexCount && DifferenceValue > ColorArray[i].TotalDifference)) {
											BestIndex = i;
											MaxSize = ColorArray[i].IndexCount;
											DifferenceValue = ColorArray[i].TotalDifference;
										}
									}

									if (BestIndex < 0) {
										break;
									}

									if (MaxSize == 1) {
										for (unsigned int i = BestIndex; i < ColorArray.size(); ++i) {
											if (ColorArray[i].IndexCount == 1) {
												tempIndexDataStream.push_back(ColorArray[i].B);
												tempIndexDataStream.push_back(ColorArray[i].G);
												tempIndexDataStream.push_back(ColorArray[i].R);

												tempColorDataStream[i] = CurrentIndex++;

												ColorArray[i].IndexCount = 0;
												ColorArray[i].IndexList.clear();
											}
										}
									} else {
										tempIndexDataStream.push_back(ColorArray[BestIndex].B);
										tempIndexDataStream.push_back(ColorArray[BestIndex].G);
										tempIndexDataStream.push_back(ColorArray[BestIndex].R);

										tempColorDataStream[BestIndex] = CurrentIndex;
										for (unsigned int i = 0; i < ColorArray[BestIndex].IndexList.size(); ++i) {
											unsigned int GroupedIndex = ColorArray[BestIndex].IndexList[i];
											if (tempColorDataStream[GroupedIndex] < 0) {
												tempColorDataStream[GroupedIndex] = CurrentIndex;
												for (unsigned int j = 0; j < ColorArray[GroupedIndex].IndexList.size(); ++j) {
													if (ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount > 0) {
														--ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount;
														if (THRESHOLD != 0) {
															ColorArray[ColorArray[GroupedIndex].IndexList[j]].TotalDifference -= ColorArray[ColorArray[GroupedIndex].IndexList[j]].DifferenceValueSet[GroupedIndex];
														}
													}
												}
												ColorArray[GroupedIndex].IndexCount = 0;
												ColorArray[GroupedIndex].IndexList.clear();
											}
										}
										ColorArray[BestIndex].IndexList.clear();
										ColorArray[BestIndex].IndexCount = 0;

										++CurrentIndex;
									}
								}
								tempIndexDataStream[0] = (char)(CurrentIndex - 1);
								ColorArray.clear();
								
								unsigned int tempIndexStreamSize = (unsigned int)IndexDataStream.size(), tempColorStreamSize = (unsigned int)ColorDataStream.size();
								IndexDataStream.resize(tempIndexDataStream.size() + tempIndexStreamSize);
								ColorDataStream.resize(tempColorDataStream.size() + tempColorStreamSize);
								for (unsigned int i = 0; i < tempIndexDataStream.size(); ++i) {
									IndexDataStream[tempIndexStreamSize + i] = tempIndexDataStream[i];
								}
								for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
									ColorDataStream[tempColorStreamSize + i] = tempColorDataStream[i];
								}
								tempIndexDataStream.clear();
								tempColorDataStream.clear();
							}
						}

						char RunValue = ColorDataStream[0];
						unsigned char LengthValue = 0;
						RunVector.push_back(RunValue);
						for (unsigned int i = 1; i < ColorDataStream.size(); ++i) {
							if (RunValue == ColorDataStream[i] && LengthValue < UCHAR_MAX) {
								++LengthValue;
							} else {
								LengthVector.push_back((char)LengthValue);
								LengthValue = 0;
								RunValue = ColorDataStream[i];
								RunVector.push_back(RunValue);
							}
						}
						LengthVector.push_back((char)LengthValue);
						ColorDataStream.clear();

						uint64_t IndexStreamSize = IndexDataStream.size();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&IndexStreamSize), sizeof (IndexStreamSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(IndexDataStream, compressed_tree_data_out_arg1);

						uint64_t RunVectorSize = RunVector.size();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&RunVectorSize), sizeof (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(RunVector, compressed_tree_data_out_arg1);

						uint64_t LengthVectorSize = LengthVector.size();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(LengthVector, compressed_tree_data_out_arg1);

						IndexDataStream.clear();
						RunVector.clear();
						LengthVector.clear();
						// encode Combination 1 end
					} else if (METHOD == 1) {
						// encode Combination 2 begin
						uint64_t Bit_Cost[2] = { 0, 0 };
						uint64_t Dist_Cost[2] = { 0, 0 };
						double RD_Cost[2] = { 0, 0 };

						vector<struct ColorData> ColorArray;
						vector<char> HeaderVector;
						vector<char> IndexDataStream, ColorDataStream;
						vector<char> RunVector, LengthVector;

						uint64_t ProbTable[4][257], tempProbTable[4][257];
						memset(ProbTable, 0, sizeof(ProbTable));
						memset(tempProbTable, 0, sizeof(tempProbTable));

						for (unsigned int idx = 0; idx < pointAvgColorDataVector.size(); idx += 3) {
							struct ColorData tempColorData;
							tempColorData.B = (unsigned char)pointAvgColorDataVector[idx];
							tempColorData.G = (unsigned char)pointAvgColorDataVector[idx + 1];
							tempColorData.R = (unsigned char)pointAvgColorDataVector[idx + 2];
							tempColorData.IndexCount = 1;
							ColorArray.push_back(tempColorData);

							if (idx + 3 == pointAvgColorDataVector.size() || ColorArray.size() == BLOCK_SIZE) {
								vector<char> tempOrigColorDataVector, tempReconColorDataVector1, tempReconColorDataVector2;
								tempOrigColorDataVector.resize(ColorArray.size() * 3);
								for (unsigned int i = 0; i < ColorArray.size(); ++i) {
									tempOrigColorDataVector[(i * 3)] = ColorArray[i].B;
									tempOrigColorDataVector[(i * 3) + 1] = ColorArray[i].G;
									tempOrigColorDataVector[(i * 3) + 2] = ColorArray[i].R;
								}
								tempReconColorDataVector1.resize(ColorArray.size() * 3);
								tempReconColorDataVector2.resize(ColorArray.size() * 3);

								// Index Codebook
								for (unsigned int i = 0; i < ColorArray.size() - 1; ++i) {
									for (unsigned int j = i + 1; j < ColorArray.size(); ++j) {
										int diff_b = (int)ColorArray[i].B - (int)ColorArray[j].B;
										int diff_g = (int)ColorArray[i].G - (int)ColorArray[j].G;
										int diff_r = (int)ColorArray[i].R - (int)ColorArray[j].R;

										double diff_rgb = pow((diff_r * diff_r) + (diff_g * diff_g) + (diff_b * diff_b), 0.5);
										if (diff_rgb <= THRESHOLD) {
											ColorArray[i].IndexList.push_back(j);
											ColorArray[j].IndexList.push_back(i);
											++ColorArray[i].IndexCount;
											++ColorArray[j].IndexCount;

											if (THRESHOLD != 0) {
												ColorArray[i].DifferenceValueSet[i] = 0;
												ColorArray[i].DifferenceValueSet[j] = (unsigned int)diff_rgb;
												ColorArray[j].DifferenceValueSet[j] = 0;
												ColorArray[j].DifferenceValueSet[i] = (unsigned int)diff_rgb;

												ColorArray[i].TotalDifference += (unsigned int)diff_rgb;
												ColorArray[j].TotalDifference += (unsigned int)diff_rgb;
											}
										}
									}
								}

								int BestIndex;
								unsigned char CurrentIndex = 0;
								unsigned int MaxSize, DifferenceValue;
								vector<char> tempIndexDataStream, tempColorDataStream;
								tempIndexDataStream.push_back(0);
								tempColorDataStream.resize(ColorArray.size());
								for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
									tempColorDataStream[i] = -1;
								}
								
								while (1) {
									BestIndex = -1;
									MaxSize = 0;

									for (unsigned int i = 0; i < ColorArray.size(); ++i) {
										if (MaxSize < ColorArray[i].IndexCount || (THRESHOLD != 0 && MaxSize == ColorArray[i].IndexCount && DifferenceValue > ColorArray[i].TotalDifference)) {
											BestIndex = i;
											MaxSize = ColorArray[i].IndexCount;
											DifferenceValue = ColorArray[i].TotalDifference;
										}
									}

									if (BestIndex < 0) {
										break;
									}

									if (MaxSize == 1) {
										for (unsigned int i = BestIndex; i < ColorArray.size(); ++i) {
											if (ColorArray[i].IndexCount == 1) {
												tempIndexDataStream.push_back(ColorArray[i].B);
												tempIndexDataStream.push_back(ColorArray[i].G);
												tempIndexDataStream.push_back(ColorArray[i].R);

												tempColorDataStream[i] = CurrentIndex++;

												tempReconColorDataVector1[i * 3] = ColorArray[i].B;
												tempReconColorDataVector1[(i * 3) + 1] = ColorArray[i].G;
												tempReconColorDataVector1[(i * 3) + 2] = ColorArray[i].R;

												ColorArray[i].IndexCount = 0;
												ColorArray[i].IndexList.clear();
											}
										}
									} else {
										tempIndexDataStream.push_back(ColorArray[BestIndex].B);
										tempIndexDataStream.push_back(ColorArray[BestIndex].G);
										tempIndexDataStream.push_back(ColorArray[BestIndex].R);

										tempColorDataStream[BestIndex] = CurrentIndex;

										tempReconColorDataVector1[BestIndex * 3] = ColorArray[BestIndex].B;
										tempReconColorDataVector1[(BestIndex * 3) + 1] = ColorArray[BestIndex].G;
										tempReconColorDataVector1[(BestIndex * 3) + 2] = ColorArray[BestIndex].R;

										for (unsigned int i = 0; i < ColorArray[BestIndex].IndexList.size(); ++i) {
											unsigned int GroupedIndex = ColorArray[BestIndex].IndexList[i];
											if (tempColorDataStream[GroupedIndex] < 0) {
												tempColorDataStream[GroupedIndex] = CurrentIndex;

												tempReconColorDataVector1[GroupedIndex * 3] = ColorArray[BestIndex].B;
												tempReconColorDataVector1[(GroupedIndex * 3) + 1] = ColorArray[BestIndex].G;
												tempReconColorDataVector1[(GroupedIndex * 3) + 2] = ColorArray[BestIndex].R;

												for (unsigned int j = 0; j < ColorArray[GroupedIndex].IndexList.size(); ++j) {
													if (ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount > 0) {
														--ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount;
														if (THRESHOLD != 0) {
															ColorArray[ColorArray[GroupedIndex].IndexList[j]].TotalDifference -= ColorArray[ColorArray[GroupedIndex].IndexList[j]].DifferenceValueSet[GroupedIndex];
														}
													}
												}
												ColorArray[GroupedIndex].IndexCount = 0;
												ColorArray[GroupedIndex].IndexList.clear();
											}
										}
										ColorArray[BestIndex].IndexList.clear();
										ColorArray[BestIndex].IndexCount = 0;

										++CurrentIndex;
									}
								}
								tempIndexDataStream[0] = (char)(CurrentIndex - 1);

								//Run-Length
								char RunValue[3];
								unsigned char LengthValue;
								vector<char> tempRunVector, tempLengthVector;

								tempRunVector.resize(3);
								tempReconColorDataVector2.resize(3);
								tempRunVector[0] = tempReconColorDataVector2[0] = RunValue[0] = ColorArray[0].B;
								tempRunVector[1] = tempReconColorDataVector2[1] = RunValue[1] = ColorArray[0].G;
								tempRunVector[2] = tempReconColorDataVector2[2] = RunValue[2] = ColorArray[0].R;
								LengthValue = 0;

								for (unsigned int i = 1; i < ColorArray.size(); ++i) {
									int diff_b = (int)((unsigned char)RunValue[0]) - (int)ColorArray[i].B;
									int diff_g = (int)((unsigned char)RunValue[1]) - (int)ColorArray[i].G;
									int diff_r = (int)((unsigned char)RunValue[2]) - (int)ColorArray[i].R;
									if (pow((diff_r * diff_r) + (diff_g * diff_g) + (diff_b * diff_b), 0.5) <= THRESHOLD && LengthValue < UCHAR_MAX) {
										++LengthValue;
									} else {
										tempLengthVector.push_back((char)LengthValue);
										RunValue[0] = ColorArray[i].B;
										RunValue[1] = ColorArray[i].G;
										RunValue[2] = ColorArray[i].R;
										tempRunVector.push_back(RunValue[0]);
										tempRunVector.push_back(RunValue[1]);
										tempRunVector.push_back(RunValue[2]);
										LengthValue = 0;
									}
									tempReconColorDataVector2.push_back(RunValue[0]);
									tempReconColorDataVector2.push_back(RunValue[1]);
									tempReconColorDataVector2.push_back(RunValue[2]);
								}
								tempLengthVector.push_back((char)LengthValue);

								ColorArray.clear();

								Bit_Cost[0] = sizeof(tempIndexDataStream.size());
								Bit_Cost[0] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempIndexDataStream, tempProbTable[0]);
								Bit_Cost[0] += sizeof(tempColorDataStream.size());
								Bit_Cost[0] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempColorDataStream, tempProbTable[1]);
								for (unsigned int i = 0; i < tempOrigColorDataVector.size(); ++i) {
									Dist_Cost[0] = ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector1[i]) * ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector1[i]);
								}
								RD_Cost[0] = (double)Dist_Cost[0] + (LAMDA[LAMDA_INDEX] * (double)Bit_Cost[0]);

								Bit_Cost[1] = sizeof(tempRunVector.size());
								Bit_Cost[1] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempRunVector, tempProbTable[2]);
								Bit_Cost[1] += sizeof(tempLengthVector.size());
								Bit_Cost[1] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempLengthVector, tempProbTable[3]);
								for (unsigned int i = 0; i < tempOrigColorDataVector.size(); ++i) {
									Dist_Cost[1] = ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector2[i]) * ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector2[i]);
								}
								RD_Cost[1] = (double)Dist_Cost[1] + (LAMDA[LAMDA_INDEX] * (double)Bit_Cost[1]);

								if (RD_Cost[0] < RD_Cost[1]) {
									/*unsigned int tempReconColorDataVectorSize = (unsigned int)ReconColorDataVector.size();
									ReconColorDataVector.resize(tempReconColorDataVector1.size() + tempReconColorDataVectorSize);
									for (unsigned int i = 0; i < tempReconColorDataVector1.size(); ++i) {
										ReconColorDataVector[tempReconColorDataVectorSize + i] = tempReconColorDataVector1[i];
									}*/

									unsigned int IndexStreamSize = (unsigned int)IndexDataStream.size();
									IndexDataStream.resize(tempIndexDataStream.size() + IndexStreamSize);
									for (unsigned int i = 0; i < tempIndexDataStream.size(); ++i) {
										IndexDataStream[IndexStreamSize + i] = tempIndexDataStream[i];
									}

									unsigned int ColorStreamSize = (unsigned int)ColorDataStream.size();
									ColorDataStream.resize(tempColorDataStream.size() + ColorStreamSize);
									for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
										ColorDataStream[ColorStreamSize + i] = tempColorDataStream[i];
									}

									HeaderVector.push_back(0);
									memcpy(ProbTable[0], tempProbTable[0], sizeof(ProbTable[0]));
									memcpy(ProbTable[1], tempProbTable[1], sizeof(ProbTable[1]));
									memcpy(tempProbTable[2], ProbTable[2], sizeof(ProbTable[2]));
									memcpy(tempProbTable[3], ProbTable[3], sizeof(ProbTable[3]));
								} else {
									/*unsigned int tempReconColorDataVectorSize = (unsigned int)ReconColorDataVector.size();
									ReconColorDataVector.resize(tempReconColorDataVector2.size() + tempReconColorDataVectorSize);
									for (unsigned int i = 0; i < tempReconColorDataVector2.size(); ++i) {
										ReconColorDataVector[tempReconColorDataVectorSize + i] = tempReconColorDataVector2[i];
									}*/

									unsigned int tempRunVectorSize = (unsigned int)RunVector.size();
									RunVector.resize(tempRunVector.size() + tempRunVectorSize);
									for (unsigned int i = 0; i < tempRunVector.size(); ++i) {
										RunVector[tempRunVectorSize + i] = tempRunVector[i];
									}

									unsigned int tempLengthVectorSize = (unsigned int)LengthVector.size();
									LengthVector.resize(tempLengthVector.size() + tempLengthVectorSize);
									for (unsigned int i = 0; i < tempLengthVector.size(); ++i) {
										LengthVector[tempLengthVectorSize + i] = tempLengthVector[i];
									}

									HeaderVector.push_back(1);
									memcpy(tempProbTable[0], ProbTable[0], sizeof(ProbTable[0]));
									memcpy(tempProbTable[1], ProbTable[1], sizeof(ProbTable[1]));
									memcpy(ProbTable[2], tempProbTable[2], sizeof(ProbTable[2]));
									memcpy(ProbTable[3], tempProbTable[3], sizeof(ProbTable[3]));
								}

								tempIndexDataStream.clear();
								tempColorDataStream.clear();
								tempRunVector.clear();
								tempLengthVector.clear();
								tempOrigColorDataVector.clear();
								tempReconColorDataVector1.clear();
								tempReconColorDataVector2.clear();
							}
						}

						uint64_t HeaderVectorSize, IndexDataStreamSize, ColorDataStreamSize, RunVectorSize, LengthVectorSize;
						HeaderVectorSize = HeaderVector.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&HeaderVectorSize), sizeof (HeaderVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(HeaderVector, compressed_tree_data_out_arg1);

						IndexDataStreamSize = IndexDataStream.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&IndexDataStreamSize), sizeof (IndexDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(IndexDataStream, compressed_tree_data_out_arg1);

						ColorDataStreamSize = ColorDataStream.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&ColorDataStreamSize), sizeof (ColorDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(ColorDataStream, compressed_tree_data_out_arg1);

						RunVectorSize = RunVector.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&RunVectorSize), sizeof (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(RunVector, compressed_tree_data_out_arg1);

						LengthVectorSize = LengthVector.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(LengthVector, compressed_tree_data_out_arg1);

						HeaderVector.clear();
						IndexDataStream.clear();
						ColorDataStream.clear();
						RunVector.clear();
						LengthVector.clear();
						// encode Combination 2 end
					} else if (METHOD == 2) {
						// encode Combination 3 begin
						uint64_t Bit_Cost[2] = { 0, 0 };
						uint64_t Dist_Cost[2] = { 0, 0 };
						double RD_Cost[2] = { 0, 0 };

						vector<struct ColorData> ColorArray;
						vector<unsigned int> HeaderVector, HeaderRunLengthVector;
						vector<char> IndexDataStream, ColorDataStream, RunVector1, LengthVector1;
						vector<char> RunVector2, LengthVector2;

						uint64_t ProbTable[4][257], tempProbTable[4][257];
						memset(ProbTable, 0, sizeof(ProbTable));
						memset(tempProbTable, 0, sizeof(tempProbTable));

						for (unsigned int idx = 0; idx < pointAvgColorDataVector.size(); idx += 3) {
							struct ColorData tempColorData;
							tempColorData.B = (unsigned char)pointAvgColorDataVector[idx];
							tempColorData.G = (unsigned char)pointAvgColorDataVector[idx + 1];
							tempColorData.R = (unsigned char)pointAvgColorDataVector[idx + 2];
							tempColorData.IndexCount = 1;
							ColorArray.push_back(tempColorData);

							if (idx + 3 == pointAvgColorDataVector.size() || ColorArray.size() == BLOCK_SIZE) {
								vector<char> tempOrigColorDataVector, tempReconColorDataVector1, tempReconColorDataVector2;
								tempOrigColorDataVector.resize(ColorArray.size() * 3);
								for (unsigned int i = 0; i < ColorArray.size(); ++i) {
									tempOrigColorDataVector[(i * 3)] = ColorArray[i].B;
									tempOrigColorDataVector[(i * 3) + 1] = ColorArray[i].G;
									tempOrigColorDataVector[(i * 3) + 2] = ColorArray[i].R;
								}
								tempReconColorDataVector1.resize(ColorArray.size() * 3);
								tempReconColorDataVector2.resize(ColorArray.size() * 3);

								// Index Codebook
								for (unsigned int i = 0; i < ColorArray.size() - 1; ++i) {
									for (unsigned int j = i + 1; j < ColorArray.size(); ++j) {
										int diff_b = (int)ColorArray[i].B - (int)ColorArray[j].B;
										int diff_g = (int)ColorArray[i].G - (int)ColorArray[j].G;
										int diff_r = (int)ColorArray[i].R - (int)ColorArray[j].R;

										double diff_rgb = pow((diff_r * diff_r) + (diff_g * diff_g) + (diff_b * diff_b), 0.5);
										if (diff_rgb <= THRESHOLD) {
											ColorArray[i].IndexList.push_back(j);
											ColorArray[j].IndexList.push_back(i);
											++ColorArray[i].IndexCount;
											++ColorArray[j].IndexCount;

											if (THRESHOLD != 0) {
												ColorArray[i].DifferenceValueSet[i] = 0;
												ColorArray[i].DifferenceValueSet[j] = (unsigned int)diff_rgb;
												ColorArray[j].DifferenceValueSet[j] = 0;
												ColorArray[j].DifferenceValueSet[i] = (unsigned int)diff_rgb;

												ColorArray[i].TotalDifference += (unsigned int)diff_rgb;
												ColorArray[j].TotalDifference += (unsigned int)diff_rgb;
											}
										}
									}
								}

								int BestIndex;
								unsigned char CurrentIndex = 0;
								unsigned int MaxSize, DifferenceValue;
								vector<char> tempIndexDataStream, tempColorDataStream;
								tempIndexDataStream.push_back(0);
								tempColorDataStream.resize(ColorArray.size());
								for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
									tempColorDataStream[i] = -1;
								}
								
								while (1) {
									BestIndex = -1;
									MaxSize = 0;

									for (unsigned int i = 0; i < ColorArray.size(); ++i) {
										if (MaxSize < ColorArray[i].IndexCount || (THRESHOLD != 0 && MaxSize == ColorArray[i].IndexCount && DifferenceValue > ColorArray[i].TotalDifference)) {
											BestIndex = i;
											MaxSize = ColorArray[i].IndexCount;
											DifferenceValue = ColorArray[i].TotalDifference;
										}
									}

									if (BestIndex < 0) {
										break;
									}

									if (MaxSize == 1) {
										for (unsigned int i = BestIndex; i < ColorArray.size(); ++i) {
											if (ColorArray[i].IndexCount == 1) {
												tempIndexDataStream.push_back(ColorArray[i].B);
												tempIndexDataStream.push_back(ColorArray[i].G);
												tempIndexDataStream.push_back(ColorArray[i].R);

												tempColorDataStream[i] = CurrentIndex++;

												tempReconColorDataVector1[i * 3] = ColorArray[i].B;
												tempReconColorDataVector1[(i * 3) + 1] = ColorArray[i].G;
												tempReconColorDataVector1[(i * 3) + 2] = ColorArray[i].R;

												ColorArray[i].IndexCount = 0;
												ColorArray[i].IndexList.clear();
											}
										}
									} else {
										tempIndexDataStream.push_back(ColorArray[BestIndex].B);
										tempIndexDataStream.push_back(ColorArray[BestIndex].G);
										tempIndexDataStream.push_back(ColorArray[BestIndex].R);

										tempColorDataStream[BestIndex] = CurrentIndex;

										tempReconColorDataVector1[BestIndex * 3] = ColorArray[BestIndex].B;
										tempReconColorDataVector1[(BestIndex * 3) + 1] = ColorArray[BestIndex].G;
										tempReconColorDataVector1[(BestIndex * 3) + 2] = ColorArray[BestIndex].R;

										for (unsigned int i = 0; i < ColorArray[BestIndex].IndexList.size(); ++i) {
											unsigned int GroupedIndex = ColorArray[BestIndex].IndexList[i];
											if (tempColorDataStream[GroupedIndex] < 0) {
												tempColorDataStream[GroupedIndex] = CurrentIndex;

												tempReconColorDataVector1[GroupedIndex * 3] = ColorArray[BestIndex].B;
												tempReconColorDataVector1[(GroupedIndex * 3) + 1] = ColorArray[BestIndex].G;
												tempReconColorDataVector1[(GroupedIndex * 3) + 2] = ColorArray[BestIndex].R;

												for (unsigned int j = 0; j < ColorArray[GroupedIndex].IndexList.size(); ++j) {
													if (ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount > 0) {
														--ColorArray[ColorArray[GroupedIndex].IndexList[j]].IndexCount;
														if (THRESHOLD != 0) {
															ColorArray[ColorArray[GroupedIndex].IndexList[j]].TotalDifference -= ColorArray[ColorArray[GroupedIndex].IndexList[j]].DifferenceValueSet[GroupedIndex];
														}
													}
												}
												ColorArray[GroupedIndex].IndexCount = 0;
												ColorArray[GroupedIndex].IndexList.clear();
											}
										}
										ColorArray[BestIndex].IndexList.clear();
										ColorArray[BestIndex].IndexCount = 0;

										++CurrentIndex;
									}
								}
								tempIndexDataStream[0] = (char)(CurrentIndex - 1);

								//Run-Length
								char tempRunValue[3];
								unsigned char tempLengthValue;
								vector<char> tempRunVector, tempLengthVector;

								tempRunVector.resize(3);
								tempReconColorDataVector2.resize(3);
								tempRunVector[0] = tempReconColorDataVector2[0] = tempRunValue[0] = ColorArray[0].B;
								tempRunVector[1] = tempReconColorDataVector2[1] = tempRunValue[1] = ColorArray[0].G;
								tempRunVector[2] = tempReconColorDataVector2[2] = tempRunValue[2] = ColorArray[0].R;
								tempLengthValue = 0;

								for (unsigned int i = 1; i < ColorArray.size(); ++i) {
									int diff_b = (int)((unsigned char)tempRunValue[0]) - (int)ColorArray[i].B;
									int diff_g = (int)((unsigned char)tempRunValue[1]) - (int)ColorArray[i].G;
									int diff_r = (int)((unsigned char)tempRunValue[2]) - (int)ColorArray[i].R;
									if (pow((diff_r * diff_r) + (diff_g * diff_g) + (diff_b * diff_b), 0.5) <= THRESHOLD && tempLengthValue < UCHAR_MAX) {
										++tempLengthValue;
									} else {
										tempLengthVector.push_back((char)tempLengthValue);
										tempRunValue[0] = ColorArray[i].B;
										tempRunValue[1] = ColorArray[i].G;
										tempRunValue[2] = ColorArray[i].R;
										tempRunVector.push_back(tempRunValue[0]);
										tempRunVector.push_back(tempRunValue[1]);
										tempRunVector.push_back(tempRunValue[2]);
										tempLengthValue = 0;
									}
									tempReconColorDataVector2.push_back(tempRunValue[0]);
									tempReconColorDataVector2.push_back(tempRunValue[1]);
									tempReconColorDataVector2.push_back(tempRunValue[2]);
								}
								tempLengthVector.push_back((char)tempLengthValue);

								ColorArray.clear();

								Bit_Cost[0] = sizeof(tempIndexDataStream.size());
								Bit_Cost[0] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempIndexDataStream, tempProbTable[0]);
								Bit_Cost[0] += sizeof(tempColorDataStream.size());
								Bit_Cost[0] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempColorDataStream, tempProbTable[1]);
								for (unsigned int i = 0; i < tempOrigColorDataVector.size(); ++i) {
									Dist_Cost[0] = ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector1[i]) * ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector1[i]);
								}
								RD_Cost[0] = (double)Dist_Cost[0] + (LAMDA[LAMDA_INDEX] * (double)Bit_Cost[0]);

								Bit_Cost[1] = sizeof(tempRunVector.size());
								Bit_Cost[1] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempRunVector, tempProbTable[2]);
								Bit_Cost[1] += sizeof(tempLengthVector.size());
								Bit_Cost[1] += entropy_coder_.encodeCharVectorToStreamForRDCost(tempLengthVector, tempProbTable[3]);
								for (unsigned int i = 0; i < tempOrigColorDataVector.size(); ++i) {
									Dist_Cost[1] = ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector2[i]) * ((int)tempOrigColorDataVector[i] - (int)tempReconColorDataVector2[i]);
								}
								RD_Cost[1] = (double)Dist_Cost[1] + (LAMDA[LAMDA_INDEX] * (double)Bit_Cost[1]);

								if (RD_Cost[0] < RD_Cost[1]) {
									/*unsigned int tempReconColorDataVectorSize = (unsigned int)ReconColorDataVector.size();
									ReconColorDataVector.resize(tempReconColorDataVector1.size() + tempReconColorDataVectorSize);
									for (unsigned int i = 0; i < tempReconColorDataVector1.size(); ++i) {
										ReconColorDataVector[tempReconColorDataVectorSize + i] = tempReconColorDataVector1[i];
									}*/

									unsigned int IndexStreamSize = (unsigned int)IndexDataStream.size();
									IndexDataStream.resize(tempIndexDataStream.size() + IndexStreamSize);
									for (unsigned int i = 0; i < tempIndexDataStream.size(); ++i) {
										IndexDataStream[IndexStreamSize + i] = tempIndexDataStream[i];
									}

									unsigned int ColorStreamSize = (unsigned int)ColorDataStream.size();
									ColorDataStream.resize(tempColorDataStream.size() + ColorStreamSize);
									for (unsigned int i = 0; i < tempColorDataStream.size(); ++i) {
										ColorDataStream[ColorStreamSize + i] = tempColorDataStream[i];
									}

									HeaderVector.push_back(0);
									memcpy(ProbTable[0], tempProbTable[0], sizeof(ProbTable[0]));
									memcpy(ProbTable[1], tempProbTable[1], sizeof(ProbTable[1]));
									memcpy(tempProbTable[2], ProbTable[2], sizeof(ProbTable[2]));
									memcpy(tempProbTable[3], ProbTable[3], sizeof(ProbTable[3]));
								} else {
									/*unsigned int tempReconColorDataVectorSize = (unsigned int)ReconColorDataVector.size();
									ReconColorDataVector.resize(tempReconColorDataVector2.size() + tempReconColorDataVectorSize);
									for (unsigned int i = 0; i < tempReconColorDataVector2.size(); ++i) {
										ReconColorDataVector[tempReconColorDataVectorSize + i] = tempReconColorDataVector2[i];
									}*/

									unsigned int tempRunVectorSize = (unsigned int)RunVector2.size();
									RunVector2.resize(tempRunVector.size() + tempRunVectorSize);
									for (unsigned int i = 0; i < tempRunVector.size(); ++i) {
										RunVector2[tempRunVectorSize + i] = tempRunVector[i];
									}

									unsigned int tempLengthVectorSize = (unsigned int)LengthVector2.size();
									LengthVector2.resize(tempLengthVector.size() + tempLengthVectorSize);
									for (unsigned int i = 0; i < tempLengthVector.size(); ++i) {
										LengthVector2[tempLengthVectorSize + i] = tempLengthVector[i];
									}

									HeaderVector.push_back(1);
									memcpy(tempProbTable[0], ProbTable[0], sizeof(ProbTable[0]));
									memcpy(tempProbTable[1], ProbTable[1], sizeof(ProbTable[1]));
									memcpy(ProbTable[2], tempProbTable[2], sizeof(ProbTable[2]));
									memcpy(ProbTable[3], tempProbTable[3], sizeof(ProbTable[3]));
								}

								tempIndexDataStream.clear();
								tempColorDataStream.clear();
								tempRunVector.clear();
								tempLengthVector.clear();
								tempOrigColorDataVector.clear();
								tempReconColorDataVector1.clear();
								tempReconColorDataVector2.clear();
							}
						}

						if (HeaderVector.size() != 0) {
							unsigned int RunValue = HeaderVector[0];
							unsigned int LengthValue = 0;
							HeaderRunLengthVector.push_back(RunValue);
							for (unsigned int i = 1; i < HeaderVector.size(); ++i) {
								if (RunValue == HeaderVector[i] && LengthValue < UINT_MAX) {
									++LengthValue;
								} else {
									HeaderRunLengthVector.push_back(LengthValue);
									RunValue = HeaderVector[i];
									LengthValue = 0;
								}
							}
							HeaderRunLengthVector.push_back(LengthValue);
						}

						if (ColorDataStream.size() != 0) {
							char RunValue = ColorDataStream[0];
							unsigned char LengthValue = 0;
							RunVector1.push_back(RunValue);
							for (unsigned int i = 1; i < ColorDataStream.size(); ++i) {
								if (RunValue == ColorDataStream[i] && LengthValue < UCHAR_MAX) {
									++LengthValue;
								} else {
									LengthVector1.push_back((char)LengthValue);
									LengthValue = 0;
									RunValue = ColorDataStream[i];
									RunVector1.push_back(RunValue);
								}
							}
							LengthVector1.push_back((char)LengthValue);
						}

						uint64_t HeaderVectorSize, IndexDataStreamSize, ColorDataStreamSize, RunVectorSize, LengthVectorSize;
						HeaderVectorSize = HeaderRunLengthVector.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&HeaderVectorSize), sizeof (HeaderVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeIntVectorToStream(HeaderRunLengthVector, compressed_tree_data_out_arg1);

						IndexDataStreamSize = IndexDataStream.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&IndexDataStreamSize), sizeof (IndexDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(IndexDataStream, compressed_tree_data_out_arg1);

						RunVectorSize = RunVector1.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&RunVectorSize), sizeof (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(RunVector1, compressed_tree_data_out_arg1);

						LengthVectorSize = LengthVector1.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(LengthVector1, compressed_tree_data_out_arg1);

						RunVectorSize = RunVector2.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&RunVectorSize), sizeof (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(RunVector2, compressed_tree_data_out_arg1);

						LengthVectorSize = LengthVector2.size ();
						compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(LengthVector2, compressed_tree_data_out_arg1);

						HeaderVector.clear();
						HeaderRunLengthVector.clear();
						IndexDataStream.clear();
						ColorDataStream.clear();
						RunVector1.clear();
						LengthVector1.clear();
						RunVector2.clear();
						LengthVector2.clear();
						// encode Combination 3 end
					}
				} else {
					point_avg_color_data_vector_size = pointAvgColorDataVector.size ();
					compressed_tree_data_out_arg1.write (reinterpret_cast<const char*> (&point_avg_color_data_vector_size), sizeof (point_avg_color_data_vector_size));
					compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream(pointAvgColorDataVector, compressed_tree_data_out_arg1);
				}
			}

			// store the number of bytes used for the color stream
			compression_performance_metrics[2] = compressed_color_data_len_ ;

			// flush output stream
			compressed_tree_data_out_arg1.flush ();

			if (!do_voxel_grid_enDecoding_)
			{
				uint64_t pointCountDataVector_size;
				uint64_t point_diff_data_vector_size;
				uint64_t point_diff_color_data_vector_size;

				// encode amount of points per voxel
				pointCountDataVector_size = point_count_data_vector_.size();
				compressed_tree_data_out_arg2.write (reinterpret_cast<const char*> (&pointCountDataVector_size), 
					sizeof (pointCountDataVector_size));
				compressed_point_data_len_ += entropy_coder_.encodeIntVectorToStream(point_count_data_vector_,
					compressed_tree_data_out_arg2);

				// encode differential point information
				std::vector<char>& point_diff_data_vector = point_coder_.getDifferentialDataVector();
				point_diff_data_vector_size = point_diff_data_vector.size ();
				compressed_tree_data_out_arg2.write (reinterpret_cast<const char*> (&point_diff_data_vector_size), sizeof (point_diff_data_vector_size));
				compressed_point_data_len_ += entropy_coder_.encodeCharVectorToStream(point_diff_data_vector, compressed_tree_data_out_arg2);

				if (cloud_with_color_)
				{
					// encode differential color information
					std::vector<char>& point_diff_color_data_vector = color_coding_type_ ?  jp_color_coder_.getDifferentialDataVector () : color_coder_.getDifferentialDataVector ();
					point_diff_color_data_vector_size = point_diff_color_data_vector.size ();
					compressed_tree_data_out_arg2.write (reinterpret_cast<const char*> (&point_diff_color_data_vector_size),
						sizeof (point_diff_color_data_vector_size));
					compressed_color_data_len_ += entropy_coder_.encodeCharVectorToStream (point_diff_color_data_vector,
						compressed_tree_data_out_arg2);
				}
			}
			// flush output stream
			compressed_tree_data_out_arg2.flush ();
		};

		/** \brief Entropy decoding of input binary stream and output to information vectors, added bitstream scalability and centroid encoding compared to  V1
		* \param compressed_tree_data_in_arg1 binary voxel input stream
		* \param compressed_tree_data_in_arg2 binary color enhancement input stream
		*/
		template<typename PointT, typename LeafT, typename BranchT, typename OctreeT> void
			OctreePointCloudCodecV2<PointT, LeafT, BranchT, OctreeT>::entropyDecoding (std::istream& compressed_tree_data_in_arg1, std::istream& compressed_tree_data_in_arg2)
		{
			uint64_t binary_tree_data_vector_size;
			uint64_t point_avg_color_data_vector_size;

			compressed_point_data_len_ = 0;
			compressed_color_data_len_ = 0;

			// decode binary octree structure
			compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&binary_tree_data_vector_size), sizeof (binary_tree_data_vector_size));
			binary_tree_data_vector_.resize(static_cast<std::size_t> (binary_tree_data_vector_size));
			compressed_point_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, binary_tree_data_vector_);

			//! new option for encoding centroids
			if(do_voxel_centroid_enDecoding_)
			{
				uint32_t l_count;

				// decode differential point information
				std::vector<char>& pointDiffDataVector = centroid_coder_.getDifferentialDataVector ();
				compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&l_count), sizeof (uint32_t));
				pointDiffDataVector.resize (static_cast<std::size_t> (l_count));
				compressed_point_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, pointDiffDataVector);
			}

			if (data_with_color_)
			{
				// decode averaged voxel color information
				std::vector<char>& point_avg_color_data_vector = color_coding_type_ ? jp_color_coder_.getAverageDataVector() : color_coder_.getAverageDataVector();

				if (PROPOSED) {
					if (METHOD == 0) {
						// decode Combination 1 begin
						uint64_t IndexStreamSize, RunVectorSize, LengthVectorSize;
						vector<char> IndexDataStream;
						vector<unsigned char> ColorDataStream;
						vector<char> RunVector, LengthVector;

						compressed_tree_data_in_arg1.read(reinterpret_cast<char*> (&IndexStreamSize), sizeof (IndexStreamSize));
						IndexDataStream.resize(static_cast<std::size_t> (IndexStreamSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, IndexDataStream);

						compressed_tree_data_in_arg1.read(reinterpret_cast<char*> (&RunVectorSize), sizeof (RunVectorSize));
						RunVector.resize(static_cast<std::size_t> (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, RunVector);

						compressed_tree_data_in_arg1.read(reinterpret_cast<char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						LengthVector.resize(static_cast<std::size_t> (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, LengthVector);

						char RunValue;
						unsigned int LengthValue;
						for (unsigned int i = 0; i < RunVector.size(); ++i) {
							RunValue = RunVector[i];
							LengthValue = (unsigned int)((unsigned char)LengthVector[i] + 1);

							unsigned int tempSize = (unsigned int)ColorDataStream.size();
							ColorDataStream.resize(static_cast<std::size_t> (tempSize + LengthValue));
							for (unsigned int j = tempSize; j < ColorDataStream.size(); ++j) {
								ColorDataStream[j] = RunValue;
							}
						}
						point_avg_color_data_vector.resize(ColorDataStream.size() * 3);
						RunVector.clear();
						LengthVector.clear();

						unsigned int IndexPosition = 0, BlockPosition = 0, DataOffset = 0;
						while (IndexPosition != IndexDataStream.size()) {
							unsigned int IndexSize = (unsigned int)((unsigned char)IndexDataStream[IndexPosition]) + 1;
							for (unsigned int i = 0; i < BLOCK_SIZE && DataOffset < point_avg_color_data_vector.size(); ++i) {
								point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned int)ColorDataStream[(BlockPosition * BLOCK_SIZE) + i] * 3) + IndexPosition + 1 ];
								point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned int)ColorDataStream[(BlockPosition * BLOCK_SIZE) + i] * 3) + IndexPosition + 2 ];
								point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned int)ColorDataStream[(BlockPosition * BLOCK_SIZE) + i] * 3) + IndexPosition + 3 ];
							}
							++BlockPosition;

							IndexPosition += (IndexSize * 3) + 1;
						}

						IndexDataStream.clear();
						ColorDataStream.clear();
						// decode Combination 1 end
					} else if (METHOD == 1) {
						// decode Combination 2 begin
						uint64_t HeaderVectorSize, IndexDataStreamSize, ColorDataStreamSize, RunVectorSize, LengthVectorSize;
						vector<char> HeaderVector;
						vector<char> IndexDataStream, ColorDataStream;
						vector<char> RunVector, LengthVector;

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&HeaderVectorSize), sizeof (HeaderVectorSize));
						HeaderVector.resize(static_cast<std::size_t> (HeaderVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, HeaderVector);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&IndexDataStreamSize), sizeof (IndexDataStreamSize));
						IndexDataStream.resize(static_cast<std::size_t> (IndexDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, IndexDataStream);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&ColorDataStreamSize), sizeof (ColorDataStreamSize));
						ColorDataStream.resize(static_cast<std::size_t> (ColorDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, ColorDataStream);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&RunVectorSize), sizeof (RunVectorSize));
						RunVector.resize(static_cast<std::size_t> (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, RunVector);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						LengthVector.resize(static_cast<std::size_t> (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, LengthVector);

						unsigned int DataSize = 0;
						DataSize += ColorDataStream.size();
						for (unsigned int i = 0; i < LengthVector.size(); ++i) {
							DataSize += (unsigned int)((unsigned char)LengthVector[i] + 1);
						}
						point_avg_color_data_vector.resize(DataSize * 3);

						unsigned int HeaderIndex = 0, DataOffset = 0;
						unsigned int IndexStreamOffset = 0, ColorStreamOffset = 0, RunVectorOffset = 0, LengthVectorOffset = 0;
						unsigned int IndexSize, BlockOffset = 0;
						while (HeaderIndex < HeaderVector.size()) {
							switch (HeaderVector[HeaderIndex++]) {
							case 0:
								IndexSize = (unsigned int)((unsigned char)IndexDataStream[IndexStreamOffset++]) + 1;
								for (unsigned int i = 0; i < BLOCK_SIZE && DataOffset < point_avg_color_data_vector.size(); ++i) {
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset];
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset + 1 ];
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset + 2 ];
								}
								IndexStreamOffset += IndexSize * 3;
								ColorStreamOffset += BLOCK_SIZE;
								break;
							case 1:
								for (unsigned int ColorDataSize = 0; ColorDataSize < BLOCK_SIZE && DataOffset < point_avg_color_data_vector.size();) {
									for (unsigned int i = 0; i <= (unsigned int)((unsigned char)LengthVector[LengthVectorOffset]); ++i) {
										point_avg_color_data_vector[DataOffset++] = RunVector[RunVectorOffset];
										point_avg_color_data_vector[DataOffset++] = RunVector[RunVectorOffset + 1];
										point_avg_color_data_vector[DataOffset++] = RunVector[RunVectorOffset + 2];
									}
									ColorDataSize += (unsigned int)((unsigned char)LengthVector[LengthVectorOffset++]) + 1;
									RunVectorOffset += 3;
								}
								break;
							default:
								printf("entropyDecoding() - Color Data Error!\n");
								exit(0);
							}
						}

						/*if (ReconColorDataVector.size() != point_avg_color_data_vector.size()) {
							printf("Error 1!\n");
							exit(0);
						}
						for (unsigned int i = 0; i < point_avg_color_data_vector.size(); ++i) {
							if (ReconColorDataVector[i] != point_avg_color_data_vector[i]) {
								printf("Error 2!\n");
								exit(0);
							}
						}*/

						HeaderVector.clear();
						IndexDataStream.clear();
						ColorDataStream.clear();
						RunVector.clear();
						LengthVector.clear();
						//ReconColorDataVector.clear();
						// decode Combination 2 end
					} else if (METHOD == 2) {
						// decode Combination 3 begin
						uint64_t HeaderVectorSize, IndexDataStreamSize, RunVectorSize, LengthVectorSize;
						vector<unsigned int> HeaderVector, HeaderRunLengthVector;
						vector<char> IndexDataStream, ColorDataStream, RunVector1, LengthVector1;
						vector<char> RunVector2, LengthVector2;

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&HeaderVectorSize), sizeof (HeaderVectorSize));
						HeaderRunLengthVector.resize(static_cast<std::size_t> (HeaderVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToIntVector(compressed_tree_data_in_arg1, HeaderRunLengthVector);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&IndexDataStreamSize), sizeof (IndexDataStreamSize));
						IndexDataStream.resize(static_cast<std::size_t> (IndexDataStreamSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, IndexDataStream);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&RunVectorSize), sizeof (RunVectorSize));
						RunVector1.resize(static_cast<std::size_t> (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, RunVector1);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						LengthVector1.resize(static_cast<std::size_t> (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, LengthVector1);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&RunVectorSize), sizeof (RunVectorSize));
						RunVector2.resize(static_cast<std::size_t> (RunVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, RunVector2);

						compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&LengthVectorSize), sizeof (LengthVectorSize));
						LengthVector2.resize(static_cast<std::size_t> (LengthVectorSize));
						compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1, LengthVector2);
						
						unsigned int HeaderValue = HeaderRunLengthVector[0];
						for (unsigned int i = 1; i < HeaderRunLengthVector.size(); ++i) {
							unsigned int HeaderVectorStart = HeaderVector.size();
							HeaderVector.resize(HeaderVector.size() + HeaderRunLengthVector[i] + 1);
							for (unsigned int j = HeaderVectorStart; j < HeaderVector.size(); ++j) {
								HeaderVector[j] = HeaderValue;
							}

							switch (HeaderValue) {
							case 0:
								HeaderValue = 1;
								break;
							case 1:
								HeaderValue = 0;
								break;
							default:
								printf("entropyDecoding() - Color Data Error!\n");
								exit(0);
							}
						}

						char RunValue;
						unsigned int LengthValue;
						for (unsigned int i = 0; i < RunVector1.size(); ++i) {
							RunValue = RunVector1[i];
							LengthValue = (unsigned int)((unsigned char)LengthVector1[i] + 1);

							unsigned int tempSize = (unsigned int)ColorDataStream.size();
							ColorDataStream.resize(static_cast<std::size_t> (tempSize + LengthValue));
							for (unsigned int j = tempSize; j < ColorDataStream.size(); ++j) {
								ColorDataStream[j] = RunValue;
							}
						}

						unsigned int DataSize = 0;
						DataSize += ColorDataStream.size();
						for (unsigned int i = 0; i < LengthVector2.size(); ++i) {
							DataSize += (unsigned int)((unsigned char)LengthVector2[i] + 1);
						}
						point_avg_color_data_vector.resize(DataSize * 3);

						unsigned int HeaderIndex = 0, DataOffset = 0;
						unsigned int IndexStreamOffset = 0, ColorStreamOffset = 0, RunVectorOffset = 0, LengthVectorOffset = 0;
						unsigned int IndexSize, BlockOffset = 0;
						while (HeaderIndex < HeaderVector.size()) {
							switch (HeaderVector[HeaderIndex++]) {
							case 0:
								IndexSize = (unsigned int)((unsigned char)IndexDataStream[IndexStreamOffset++]) + 1;
								for (unsigned int i = 0; i < BLOCK_SIZE && DataOffset < point_avg_color_data_vector.size(); ++i) {
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset];
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset + 1 ];
									point_avg_color_data_vector[DataOffset++] = IndexDataStream[ ((unsigned char)ColorDataStream[ColorStreamOffset + i] * 3) + IndexStreamOffset + 2 ];
								}
								IndexStreamOffset += IndexSize * 3;
								ColorStreamOffset += BLOCK_SIZE;
								break;
							case 1:
								for (unsigned int ColorDataSize = 0; ColorDataSize < BLOCK_SIZE && DataOffset < point_avg_color_data_vector.size();) {
									for (unsigned int i = 0; i <= (unsigned int)((unsigned char)LengthVector2[LengthVectorOffset]); ++i) {
										point_avg_color_data_vector[DataOffset++] = RunVector2[RunVectorOffset];
										point_avg_color_data_vector[DataOffset++] = RunVector2[RunVectorOffset + 1];
										point_avg_color_data_vector[DataOffset++] = RunVector2[RunVectorOffset + 2];
									}
									ColorDataSize += (unsigned int)((unsigned char)LengthVector2[LengthVectorOffset++]) + 1;
									RunVectorOffset += 3;
								}
								break;
							default:
								printf("entropyDecoding() - Color Data Error!\n");
								exit(0);
							}
						}

						/*if (ReconColorDataVector.size() != point_avg_color_data_vector.size()) {
							printf("Error 1!\n");
							exit(0);
						}
						for (unsigned int i = 0; i < point_avg_color_data_vector.size(); ++i) {
							if (ReconColorDataVector[i] != point_avg_color_data_vector[i]) {
								printf("Error 2!\n");
								exit(0);
							}
						}*/

						HeaderVector.clear();
						HeaderRunLengthVector.clear();
						IndexDataStream.clear();
						ColorDataStream.clear();
						RunVector1.clear();
						LengthVector1.clear();
						RunVector2.clear();
						LengthVector2.clear();
						//ReconColorDataVector.clear();
						// decode Combination 3 end
					}
					compressed_tree_data_in_arg1.read (reinterpret_cast<char*> (&point_avg_color_data_vector_size), sizeof (point_avg_color_data_vector_size));
					point_avg_color_data_vector.resize (static_cast<std::size_t> (point_avg_color_data_vector_size));
					compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector(compressed_tree_data_in_arg1,	point_avg_color_data_vector);
				}
			}

			//! check if the enhancement layer has been received
			compressed_tree_data_in_arg2.peek();
			if(compressed_tree_data_in_arg2.good() && (! compressed_tree_data_in_arg2.eof()) )
				do_voxel_grid_enDecoding_ = false;
			else
				do_voxel_grid_enDecoding_ = true;

			if (!do_voxel_grid_enDecoding_)
			{
				uint64_t point_count_data_vector_size;
				uint64_t point_diff_data_vector_size;
				uint64_t point_diff_color_data_vector_size;

				// decode amount of points per voxel
				compressed_tree_data_in_arg2.read (reinterpret_cast<char*> (&point_count_data_vector_size), sizeof (point_count_data_vector_size));
				point_count_data_vector_.resize (static_cast<std::size_t> (point_count_data_vector_size));
				compressed_point_data_len_ += entropy_coder_.decodeStreamToIntVector(compressed_tree_data_in_arg2, point_count_data_vector_);
				point_count_data_vector_iterator_ = point_count_data_vector_.begin ();

				// decode differential point information
				std::vector<char>& pointDiffDataVector = point_coder_.getDifferentialDataVector();
				compressed_tree_data_in_arg2.read (reinterpret_cast<char*> (&point_diff_data_vector_size), sizeof (point_diff_data_vector_size));
				pointDiffDataVector.resize (static_cast<std::size_t> (point_diff_data_vector_size));
				compressed_point_data_len_ += entropy_coder_.decodeStreamToCharVector (compressed_tree_data_in_arg2, pointDiffDataVector);
				if (data_with_color_)
				{
					// decode differential color information
					std::vector<char>& pointDiffColorDataVector = color_coder_.getDifferentialDataVector ();
					compressed_tree_data_in_arg2.read (reinterpret_cast<char*> (&point_diff_color_data_vector_size), sizeof (point_diff_color_data_vector_size));
					pointDiffColorDataVector.resize (static_cast<std::size_t> (point_diff_color_data_vector_size));
					compressed_color_data_len_ += entropy_coder_.decodeStreamToCharVector (compressed_tree_data_in_arg2, pointDiffColorDataVector);
				}
			}
		};
	}
}
#endif
