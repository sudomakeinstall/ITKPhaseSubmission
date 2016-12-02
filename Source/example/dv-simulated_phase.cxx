/////////////////////////
// Contributed Classes //
/////////////////////////

#include "itkWrapPhaseFunctor.h"
#include "itkPhaseExamplesImageSource.h"
#include "itkItohPhaseUnwrappingImageFilter.h"
#include "itkPhaseResidueImageFilter.h"
#include "itkHelmholtzDecompositionImageFilter.h"
#include "itkQualityGuidedPhaseUnwrappingImageFilter.h"
#include "itkDCTPhaseUnwrappingImageFilter.h"
#include "itkPhaseQualityImageFilter.h"

///////////////////
// Other Classes //
///////////////////

#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileWriter.h"
#include "vnl/vnl_math.h"
#include "itkLineIterator.h"
#include <fstream>
#include "itkTimeProbe.h"

//////////////
// Typedefs //
//////////////

// General
const unsigned int Dimension = 2;
typedef double                                   WorkPixelType;
typedef unsigned short                           OutputPixelType;
typedef itk::Image< WorkPixelType, Dimension >   WorkImageType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
typedef itk::LineIterator< WorkImageType >       ItType;

// ITK Filters
typedef itk::ImageFileReader< WorkImageType >                              ReaderType;
typedef itk::RescaleIntensityImageFilter< WorkImageType, OutputImageType > RescaleType;
typedef itk::RescaleIntensityImageFilter< WorkImageType >                  RescalePIType;
typedef itk::ImageFileWriter< OutputImageType >                            WriterType;

// Contributed Filters
typedef itk::PhaseExamplesImageSource< WorkImageType >                ExampleType;
typedef itk::ItohPhaseUnwrappingImageFilter< WorkImageType >          ItohType;
typedef itk::PhaseResidueImageFilter< WorkImageType >                 ResidueType;
typedef itk::HelmholtzDecompositionImageFilter< WorkImageType >       HDType;
typedef itk::PhaseQualityImageFilter< WorkImageType >                 QualityType;
typedef itk::QualityGuidedPhaseUnwrappingImageFilter< WorkImageType > QGUnwrapType;
typedef itk::DCTPhaseUnwrappingImageFilter< WorkImageType >           DCTUnwrapType;

int main( int argc, char * argv[] )
{

  /////////////////////////////////////////
  // Introduction to the ITKPhase Module //
  /////////////////////////////////////////

  itk::Functor::WrapPhaseFunctor< double > wrapFunc;
  
  std::cout << wrapFunc( 3.0 ) << std::endl; // 3
  std::cout << wrapFunc( 0.0 ) << std::endl; // 0
  std::cout << wrapFunc( -3.0 ) << std::endl; // -3
  std::cout << wrapFunc( 1.0 + vnl_math::pi ) << std::endl; // -2.14159
  std::cout << wrapFunc( -vnl_math::pi - 1.0 ) << std::endl; // 2.14159
  
  ExampleType::Pointer phase = ExampleType::New();
  RescaleType::Pointer rescale = RescaleType::New();
  WriterType::Pointer writer = WriterType::New();
	
  // The default simulated phase example is a wrapped phase ramp.
  phase->Update();
  rescale->SetInput( phase->GetOutput() ); // Rescaled for visualization purposes only
  writer->SetInput( rescale->GetOutput() );
  writer->SetFileName( "../data/2/00a_ramp_wrapped.png" );
  writer->Update();
  
  // The SetWrap(bool) method determines whether the image is wrapped or unwrapped
  phase->SetWrap( false );
  phase->Update();
  writer->SetFileName( "../data/2/00b_ramp_unwrapped.png" );
  writer->Update();
	
  // Add a patch of noise with the SetNoise(bool) method
  // For reproducibility, set the seed using SetSeed(unsigned int)
  phase->SetWrap( true );
  phase->SetNoise( true );
  phase->SetNoiseSeed( 0 );
  phase->Update();
  writer->SetFileName( "../data/2/00c_noise_wrapped.png" );
  writer->Update();
	
  // Get the unwrapped, noisy example
  phase->SetWrap( false );
  phase->Update();
  writer->SetFileName( "../data/2/00d_noise_unwrapped.png" );
  writer->Update();
	
  //////////////////////////////////////
  // Phase Quality and Path Following //
  /////////////////////////////////////

  // The Itoh phase unwrapping algorithm is the simplest possible path-following algorithm
  // This algorithm proceeds in one direction only
  ItohType::Pointer itoh1 = ItohType::New();
  ItohType::Pointer itoh2 = ItohType::New();

  // Set up the example phase image to be wrapped and without a noise patch.
  phase->SetWrap( true );
  phase->SetNoise( true );
  
  itoh1->SetDirection( 0 ); // first x
  itoh2->SetDirection( 1 ); // then y

  itoh1->SetInput( phase->GetOutput() );
  itoh2->SetInput( itoh1->GetOutput() );

  rescale->SetInput( itoh2->GetOutput() );
  writer->SetFileName( "../data/3/01a_noise_xy.png" );
  writer->Update();
	
  itoh1->SetDirection( 1 );
  itoh2->SetDirection( 0 );

  writer->SetFileName( "../data/3/01b_noise_yx.png" );
  writer->Update();

  // Residues
  ResidueType::Pointer residue = ResidueType::New();
  QualityType::Pointer qual = QualityType::New();
  HDType::Pointer hd = HDType::New();
	
  // Residue
  residue->SetInput( phase->GetOutput() );
  rescale->SetInput( residue->GetOutput() );
  writer->SetFileName( "../data/3/02a_residue.png" );
  writer->Update();
  
  // Helmholtz Decomposition
  hd->SetInput( phase->GetOutput() );
  
  rescale->SetInput( hd->GetIrrotational() );
  writer->SetFileName( "../data/3/02b_irrotational.png" );
  writer->Update();
  
  rescale->SetInput( hd->GetRotational() );
  writer->SetFileName( "../data/3/02c_rotational.png" );
  writer->Update();
  
  // Phase derivative Variance
  qual->SetInput( phase->GetOutput() );
  rescale->SetInput( qual->GetOutput() );
  writer->SetFileName( "../data/3/02d_quality.png" );
  writer->Update();

  ///////////////////////////////
  // Filters for HARP/SWI Data //
  ///////////////////////////////

  ReaderType::Pointer reader = ReaderType::New();
  RescalePIType::Pointer rescalePI = RescalePIType::New();
  QGUnwrapType::Pointer qgUnwrap = QGUnwrapType::New();
  DCTUnwrapType::Pointer dctUnwrap = DCTUnwrapType::New();
  WorkImageType::IndexType truePhase;
  
  ///////////////
  // HARP data //
  ///////////////

//  std::string input = "wrapped_swi.png";
  std::string input = "wrapped_harp.png";
  std::string dataDir = "../data/";
  std::string inputDir = dataDir + "input/";
//  std::string outputDir = dataDir + "swi/"; 
  std::string outputDir = dataDir + "harp/";

  reader->SetFileName( inputDir + input );
  
  rescalePI->SetInput( reader->GetOutput() );
  rescalePI->SetOutputMaximum( vnl_math::pi );
  rescalePI->SetOutputMinimum( -1*vnl_math::pi );

  // Residue
  residue->SetInput( rescalePI->GetOutput() );
  rescale->SetInput( residue->GetOutput() );
  writer->SetInput( rescale->GetOutput() );
  writer->SetFileName( outputDir + "03a_residue.png" );
  writer->Update();
	
  // HD
  hd->SetInput( rescalePI->GetOutput() );
  rescale->SetInput( hd->GetIrrotational() );
  writer->SetFileName( outputDir + "03b_irrotational.png" );
  writer->Update();
  
  rescale->SetInput( hd->GetRotational() );
  writer->SetFileName( outputDir + "03c_rotational.png" );
  writer->Update();

  // Phase Derivative Variance
  qual->SetInput( rescalePI->GetOutput() );
  rescale->SetInput( qual->GetOutput() );
  writer->SetFileName( outputDir + "03d_quality.png" );
  writer->Update();

  // Quality Guided Phase
  truePhase[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] / 2;
  truePhase[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] / 2;

  qgUnwrap->SetInput( rescalePI->GetOutput() );
  qgUnwrap->SetTruePhase( truePhase );

  itk::TimeProbe clock;
  clock.Start();
  qgUnwrap->Update();
  clock.Stop();
  std::cout << "Quality-Guided Time Elapsed: " << clock.GetTotal() << std::endl;
  rescale->SetInput( qgUnwrap->GetPhase() );
  writer->SetFileName( outputDir + "03e_qg_unwrap.png" );
  writer->Update();
	
  dctUnwrap->SetInput( rescalePI->GetOutput() );
  clock.Reset();
  clock.Start();
  dctUnwrap->Update();
  clock.Stop();
  std::cout << "DCT Time Elapsed: " << clock.GetTotal() << std::endl;
  rescale->SetInput( dctUnwrap->GetOutput() );
  writer->SetFileName( outputDir + "03f_dct_unwrap.png" );
  writer->Update();
  
  std::ofstream file;
  file.open( (outputDir + "congruence.csv").c_str() );
  file << "Index,Wrapped,Quality,DCT\n";
 
  WorkImageType::IndexType start;
  start[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0]/2;
  start[1] = 0;

  WorkImageType::IndexType end;
  end[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0]/2;
  end[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1] - 1;

  ItType wrappedIt( rescalePI->GetOutput(), start, end );
  ItType qgIt( qgUnwrap->GetPhase(), start, end );
  ItType dctIt( dctUnwrap->GetOutput(), start, end );

  wrappedIt.GoToBegin();
  qgIt.GoToBegin();
  dctIt.GoToBegin();

  while (!wrappedIt.IsAtEnd()) {
    
    file << wrappedIt.GetIndex()[1] << ",";
    file << wrappedIt.Get() << ",";
    file << qgIt.Get() << ",";
    file << dctIt.Get() << "\n";

    ++wrappedIt;
    ++qgIt;
    ++dctIt;

  }

  file.close();

  return EXIT_SUCCESS;
 
}
