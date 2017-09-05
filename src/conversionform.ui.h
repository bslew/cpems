/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/




void ConversionForm::init()
{
    //        numberLineEdit->setValidator( new QDoubleValidator( numberLineEdit ) );
    lineEdit1->setText( "10" );
    //        convert();
    //        numberLineEdit->selectAll();
}



void ConversionForm::sum_stuff()
{
    double input1 = lineEdit1->text().toDouble();
    double input2 = lineEdit2->text().toDouble();
    simple_algebra goalgebra;
    double result =     goalgebra.get_sum(input1,input2);
    lineEdit3->setText( QString::number( result, 'f', 3 ) );
        
}











void ConversionForm::feedbackresult()
{
    double tmp = lineEdit3->text().toDouble();
     lineEdit2->setText( QString::number( tmp, 'f', 3 ) );    

}
