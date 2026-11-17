#include "bout/bout.hxx"
#include "bout/field_factory.hxx"
#include "bout/invert_pardiv.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  {
    auto inv = InvertParDiv::create(&Options::root()["pardiv"]);
    inv->setCoefA(1.0);
    inv->setCoefB(FieldFactory::get()->create3D("1 + 0.1*sin(z)"));
    Field3D f = FieldFactory::get()->create3D("sin(y) + cos(z)");
    Field3D x = inv->solve(f);
  }
  BoutFinalise();
  return 0;
}
